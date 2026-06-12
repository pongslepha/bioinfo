############################################################
## GSE39872 LIN28 CLIP-seq + RNA-seq (knockdown) pathway analysis
##
## GSE39872 (Wilbert et al., Mol Cell 2012,
## "LIN28 binds messenger RNAs at the loops of stem-loop
## structures and regulates RNA target abundance") 은 03/04
## 스크립트와 동일한 LIN28 단백질을 다룬다. 다만 human (H9 hESC,
## hg18) 데이터라서 mouse 데이터와 다음 차이가 있다.
##
## 03 (mouse, GSE) / 04 (GSE134033, mouse) 와 비교:
##   - 종(species)    : human  -> org.Hs.eg.db 사용 (mouse: org.Mm.eg.db)
##   - x축 CLIP signal : gene-level FPKM/count 가 아니라 "CLIP cluster
##                       BED"(유전자 좌표 없음)다. 따라서 hg18 knownGene
##                       좌표와 strand-specific 하게 overlap 시켜 유전자별
##                       CLIP read count 를 직접 만든다.
##   - y축 functional  : ribosome profiling 이 없으므로 04 와 동일하게
##                       "LIN28 knockdown vs control knockdown" mRNA
##                       발현 변화로 둔다. (LIN28 은 표적 mRNA 안정성을
##                       조절하므로 KD 시 표적 발현이 변한다.)
##   - localization    : mouselocalization 파일은 ENSMUSG(mouse) 기준이라
##                       human 유전자에 직접 join 되지 않는다. 따라서
##                       localization group 분석은 생략하고, 03/04 의
##                       핵심인 GO(CC/BP/MF) enrichment 분석만 수행한다.
##
## 분석 틀 자체는 03/04 와 동일:
##   x = LIN28 CLIP enrichment (log2),  y = LIN28-KD mRNA 발현 변화 (log2)
##   -> GO term 별 group-wise Wilcoxon -> bubble plot + lollipop dotplot
##
## Input:
##   00.data/GEO/GSE39872/GSM980593_LIN28ES_CLIPseq_clusters.BED
##       LIN28 CLIP cluster (H9 ES cells, hg18). BED col5 = cluster read count.
##   00.data/GEO/GSE39872/GSM980596_LIN28KD_H9_rnaseq_RPKM.txt   (case  = LIN28 KD)
##   00.data/GEO/GSE39872/GSM980597_controlKD_H9_rnaseq_RPKM.txt (control = control KD)
##       UCSC_KnownGene_ID(uc####.v) + RPKM
##   00.data/hg18_ucsc/knownGene.txt.gz        (uc id -> hg18 좌표)
##   00.data/hg18_ucsc/knownToLocusLink.txt.gz (uc id -> Entrez)
##
## 실행: lab conda env R (org.Hs.eg.db / GO.db / GenomicRanges 보유)
##   ~/miniconda3/envs/lab/bin/Rscript xx.script/05.pathway_analysis_GSE39872.R
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(ggrepel)
})

############################################################
## 0. Paths  (03 스크립트와 동일하게 스크립트 위치 기준으로 setwd)
############################################################
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}

setwd(get_script_dir())
cat("Working directory set to:", getwd(), "\n")

data_dir   <- "../00.data"
geo_dir    <- file.path(data_dir, "GEO", "GSE39872")
hg18_dir   <- file.path(data_dir, "hg18_ucsc")
output_dir <- "../01.result/own_analysis/GSE39872_analysis"

clip_file       <- file.path(geo_dir, "GSM980593_LIN28ES_CLIPseq_clusters.BED")
rna_case_file   <- file.path(geo_dir, "GSM980596_LIN28KD_H9_rnaseq_RPKM.txt")   # LIN28 KD
rna_ctrl_file   <- file.path(geo_dir, "GSM980597_controlKD_H9_rnaseq_RPKM.txt") # control KD
knownGene_file  <- file.path(hg18_dir, "knownGene.txt.gz")
uc2entrez_file  <- file.path(hg18_dir, "knownToLocusLink.txt.gz")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pc <- 1  # log 계산용 pseudocount (03/04 와 동일)

############################################################
## 1. hg18 knownGene 좌표 + uc->Entrez 매핑 로드
############################################################
## CLIP cluster(BED) 를 유전자에 배정하려면 유전자 좌표가 필요하다.
## hg18 UCSC knownGene flat 파일을 직접 사용한다 (build = hg18, RNA-seq
## uc id 와 100% 일치 확인됨).

kg <- fread(
  cmd = paste("zcat", shQuote(knownGene_file)),
  header = FALSE, sep = "\t",
  select = 1:5,
  col.names = c("uc", "chrom", "strand", "txStart", "txEnd")
)

uc2entrez <- fread(
  cmd = paste("zcat", shQuote(uc2entrez_file)),
  header = FALSE, sep = "\t",
  col.names = c("uc", "entrez")
) %>%
  mutate(entrez = as.character(entrez)) %>%
  distinct(uc, .keep_all = TRUE)

## 표준 염색체만 유지 (chr1..22, X, Y) — random/hap contig 제외
std_chrom <- paste0("chr", c(1:22, "X", "Y"))

kg <- kg %>%
  filter(chrom %in% std_chrom) %>%
  left_join(uc2entrez, by = "uc") %>%
  filter(!is.na(entrez))

## uc 별 transcript-span GRanges (BED/UCSC 0-based start -> 1-based +1).
## CLIP read 는 mRNA 와 같은 strand 에 찍히므로 strand 정보를 유지해
## strand-specific overlap 에 사용한다.
kg_gr <- GRanges(
  seqnames = kg$chrom,
  ranges   = IRanges(start = kg$txStart + 1, end = kg$txEnd),
  strand   = kg$strand
)
mcols(kg_gr)$entrez <- kg$entrez

############################################################
## 2. CLIP cluster(BED) 로드 -> 유전자별 CLIP read count
############################################################
## BED 컬럼: chrom,start,end, col4(cluster id), col5(=cluster read count),
##           strand, thickStart, thickEnd, rgb
## col5 를 CLIP 신호(읽힌 read 수)로 사용한다.

clip_raw <- fread(
  clip_file, header = FALSE, sep = "\t",
  select = 1:6,
  col.names = c("chrom", "start", "end", "cluster_id", "reads", "strand")
) %>%
  filter(chrom %in% std_chrom)

clip_gr <- GRanges(
  seqnames = clip_raw$chrom,
  ranges   = IRanges(start = clip_raw$start + 1, end = clip_raw$end),
  strand   = clip_raw$strand
)
mcols(clip_gr)$reads <- clip_raw$reads

## strand-specific overlap: cluster -> uc transcript
hits <- findOverlaps(clip_gr, kg_gr, ignore.strand = FALSE)

## (cluster, entrez) 단위로 묶어, 한 cluster 가 같은 유전자의 여러 isoform
## 에 걸쳐도 그 유전자에 read 수를 한 번만 더한다.
clip_hit_df <- tibble(
  cluster = queryHits(hits),
  entrez  = mcols(kg_gr)$entrez[subjectHits(hits)],
  reads   = mcols(clip_gr)$reads[queryHits(hits)]
) %>%
  distinct(cluster, entrez, .keep_all = TRUE)

clip_gene <- clip_hit_df %>%
  group_by(entrez) %>%
  summarise(clip_reads = sum(reads), .groups = "drop")

## 라이브러리 크기(전체 cluster read 합)로 CPM 정규화
total_clip_reads <- sum(clip_raw$reads)
clip_gene <- clip_gene %>%
  mutate(clip_cpm = clip_reads / total_clip_reads * 1e6)

message(sprintf("CLIP clusters: %d, assigned to %d genes (%.1f%% of read mass)",
                nrow(clip_raw), nrow(clip_gene),
                100 * sum(clip_gene$clip_reads) / total_clip_reads))

############################################################
## 3. RNA-seq(RPKM) 로드 -> 유전자(Entrez) 수준으로 집계
############################################################
## uc id 별 RPKM 을 Entrez 로 매핑한 뒤, 유전자별로 isoform RPKM 을 합산해
## gene-level 발현으로 사용한다. control KD 와 LIN28 KD 를 uc 로 join.

read_rpkm <- function(path, valname) {
  fread(path, header = TRUE, sep = "\t",
        col.names = c("uc", valname))
}

rna_ctrl <- read_rpkm(rna_ctrl_file, "rpkm_control")  # control KD
rna_case <- read_rpkm(rna_case_file, "rpkm_case")     # LIN28 KD

rna_df <- full_join(rna_ctrl, rna_case, by = "uc") %>%
  mutate(
    rpkm_control = replace_na(rpkm_control, 0),
    rpkm_case    = replace_na(rpkm_case, 0)
  ) %>%
  left_join(uc2entrez, by = "uc") %>%
  filter(!is.na(entrez)) %>%
  group_by(entrez) %>%
  summarise(
    rna_control = sum(rpkm_control),  # control KD gene-level RPKM
    rna_case    = sum(rpkm_case),     # LIN28 KD gene-level RPKM
    .groups = "drop"
  )

## control KD 발현도 CPM 스케일로 바꿔, CLIP_cpm 과 같은 척도에서 비율을 잡는다.
rna_df <- rna_df %>%
  mutate(rna_control_cpm = rna_control / sum(rna_control) * 1e6)

############################################################
## 4. CLIP + RNA join 및 두 지표 계산
############################################################
## clip_enrichment : LIN28 이 그 전사체에 얼마나 결합하는가 (03/04 x축 개념)
##     = log2( (CLIP_cpm + pc) / (control RNA_cpm + pc) )
##     CLIP 신호가 없으면 clip_cpm = 0 (결합 없음) 으로 둔다.
## expr_change     : LIN28 knockdown 시 mRNA 수준 변화 (04 y축과 동일)
##     = log2( (LIN28KD RPKM + pc) / (control KD RPKM + pc) )

gene_stat <- rna_df %>%
  left_join(clip_gene %>% dplyr::select(entrez, clip_cpm), by = "entrez") %>%
  mutate(
    clip_cpm        = replace_na(clip_cpm, 0),
    clip_enrichment = log2((clip_cpm + pc) / (rna_control_cpm + pc)),
    expr_change     = log2((rna_case + pc) / (rna_control + pc))
  )

############################################################
## 5. gene symbol 부착 + 기본 필터링
############################################################
## 발현이 거의 없는 유전자 제거 + 03/04 와 동일하게 histone 유전자
## (poly(A) 정규화 artifact 가능성) 제거.

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gene_stat$entrez),
  keytype = "ENTREZID",
  columns = c("SYMBOL")
) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gene_stat <- gene_stat %>%
  left_join(gene_map, by = c("entrez" = "ENTREZID")) %>%
  filter(
    (rna_control + rna_case) >= 1,
    is.finite(clip_enrichment),
    is.finite(expr_change),
    !is.na(SYMBOL),
    !str_detect(SYMBOL, "^HIST")
  ) %>%
  dplyr::rename(symbol = SYMBOL)

message(sprintf("Genes after filtering: %d (with CLIP signal: %d)",
                nrow(gene_stat), sum(gene_stat$clip_cpm > 0)))

############################################################
## 6. gene-level table 저장
############################################################
write_csv(
  gene_stat,
  file.path(output_dir, "gene_level_CLIP_RNA_table.csv")
)

############################################################
## 7. GO enrichment Figure 5A-like plots (CC / BP / MF)
############################################################
## 03/04 의 14 / 9 단계와 동일. Entrez -> GO 매핑 후 ontology 별
## group-wise Wilcoxon (그룹 유전자 vs 나머지) -> 통계표 + bubble + dotplot.

go_map_all <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gene_stat$entrez),
  keytype = "ENTREZID",
  columns = c("GO", "ONTOLOGY")
) %>%
  filter(!is.na(GO)) %>%
  distinct(ENTREZID, GO, .keep_all = TRUE) %>%
  dplyr::rename(entrez = ENTREZID)

all_genes_go <- gene_stat %>%
  dplyr::select(entrez, clip_enrichment, expr_change) %>%
  distinct()

ont_full_name <- c(
  CC = "Cellular Component",
  BP = "Biological Process",
  MF = "Molecular Function",
  ALL = "All ontologies (CC + BP + MF)"
)

## 한 ontology 에 대해 group-wise Wilcoxon 검정 -> 통계표 + plot 저장.
## ont == "ALL" 이면 세 ontology 를 구분 없이 합쳐서 분석한다. (원본 03 과 동일)
run_go_ontology <- function(ont) {
  message("== GO ", ont, " analysis ==")

  go_map <- if (ont == "ALL") {
    go_map_all
  } else {
    go_map_all %>% filter(ONTOLOGY == ont)
  }

  go_term_names <- AnnotationDbi::select(
    GO.db,
    keys = unique(go_map$GO),
    keytype = "GOID",
    columns = c("TERM")
  ) %>%
    distinct(GOID, .keep_all = TRUE)

  dat_go <- go_map %>% inner_join(all_genes_go, by = "entrez")

  test_one_go <- function(go_id) {
    genes_in <- dat_go %>% filter(GO == go_id) %>% pull(entrez) %>% unique()
    in_df  <- all_genes_go %>% filter(entrez %in% genes_in)
    out_df <- all_genes_go %>% filter(!entrez %in% genes_in)
    if (nrow(in_df) < 10 || nrow(out_df) < 10) return(NULL)

    tibble(
      GO                   = go_id,
      n_genes              = nrow(in_df),
      mean_clip_enrichment = mean(in_df$clip_enrichment, na.rm = TRUE),
      mean_expr_change     = mean(in_df$expr_change, na.rm = TRUE),
      clip_p = wilcox.test(in_df$clip_enrichment, out_df$clip_enrichment)$p.value,
      expr_p = wilcox.test(in_df$expr_change,     out_df$expr_change)$p.value
    )
  }

  go_stats <- map_dfr(unique(dat_go$GO), test_one_go) %>%
    left_join(go_term_names, by = c("GO" = "GOID")) %>%
    mutate(
      ontology     = ont,
      clip_fdr     = p.adjust(clip_p, method = "BH"),
      expr_fdr     = p.adjust(expr_p, method = "BH"),
      combined_fdr = pmin(clip_fdr, expr_fdr, na.rm = TRUE),
      neglog10_fdr = -log10(combined_fdr)
    ) %>%
    filter(n_genes >= 10, n_genes <= 500, combined_fdr < 0.05) %>%
    arrange(combined_fdr)

  write_csv(
    go_stats,
    file.path(output_dir, paste0("GO_", ont, "_enrichment_stats.csv"))
  )

  ## ----- Figure 5A-like 2D bubble plot -----
  label_go <- go_stats %>% slice_head(n = 25)

  p_bubble <- ggplot(
    go_stats %>% arrange(neglog10_fdr),
    aes(x = mean_clip_enrichment, y = mean_expr_change)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(size = n_genes, color = neglog10_fdr), alpha = 0.75) +
    geom_text_repel(
      data = label_go,
      aes(label = paste0(TERM, " (", n_genes, ")")),
      size = 3, max.overlaps = Inf
    ) +
    scale_size_continuous(range = c(2, 14), name = "Number of genes") +
    scale_color_gradient(low = "khaki1", high = "red4", name = "-log10(FDR)") +
    labs(
      title = paste0("GSE39872 GO ", ont_full_name[[ont]],
                     " analysis (LIN28 CLIP vs LIN28-KD expression)"),
      x = "Mean LIN28 CLIP enrichment (log2)",
      y = "Mean log2 fold change (LIN28 KD / control KD)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_line(color = "grey85", linetype = "dashed"),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(output_dir, paste0("GO_", ont, "_Figure5A_like.png")),
    p_bubble, width = 12, height = 7, dpi = 900
  )

  ## ----- Dot (lollipop) plot ----- (03/04 와 동일한 스타일)
  n_top <- 10

  ## facet_wrap(ncol = 2) 은 panel 을 행 우선(row-major)으로 채우므로,
  ## level 순서를 row1=Negative, row2=Positive / col1=CLIP, col2=Expression 으로 맞춘다.
  metric_levels <- c("CLIP enrichment", "Expression change (KD)")  # column 순서
  dir_levels    <- c("Negative", "Positive")                       # row 순서

  go_long <- bind_rows(
    go_stats %>% transmute(
      TERM, n_genes, metric = "CLIP enrichment",
      fdr = clip_fdr, effect = mean_clip_enrichment
    ),
    go_stats %>% transmute(
      TERM, n_genes, metric = "Expression change (KD)",
      fdr = expr_fdr, effect = mean_expr_change
    )
  ) %>%
    filter(fdr < 0.05) %>%
    mutate(
      neglog10_fdr = -log10(pmax(fdr, 1e-320)),
      direction    = if_else(effect >= 0, "Positive", "Negative")
    ) %>%
    group_by(metric, direction) %>%
    slice_min(fdr, n = n_top, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      facet_lab = factor(
        paste0(metric, "\n", direction),
        ## row(direction) 별로 column(metric) 을 채워 row-major 순서를 만든다.
        levels = as.vector(t(outer(dir_levels, metric_levels,
                                    function(d, m) paste0(m, "\n", d))))
      )
    ) %>%
    arrange(facet_lab, neglog10_fdr) %>%
    mutate(term_ord = factor(row_number(), levels = row_number(), labels = TERM))

  p_dot <- ggplot(go_long, aes(x = neglog10_fdr, y = term_ord)) +
    geom_segment(aes(x = 0, xend = neglog10_fdr, yend = term_ord, color = direction)) +
    geom_point(aes(size = n_genes, color = direction)) +
    facet_wrap(~ facet_lab, scales = "free_y", ncol = 2) +
    scale_color_manual(
      values = c(Positive = "firebrick", Negative = "steelblue"),
      name = "Effect direction"
    ) +
    scale_size_continuous(range = c(2, 9), name = "Number of genes") +
    labs(
      title = paste0(
        "GSE39872 GO ", ont_full_name[[ont]],
        " enrichment (LIN28 CLIP vs LIN28-KD expression; top ", n_top,
        " per direction)"
      ),
      x = "-log10(FDR)", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(output_dir, paste0("GO_", ont, "_dotplot.png")),
    p_dot, width = 16, height = 8, dpi = 900
  )

  go_stats
}

go_results <- lapply(c("CC", "BP", "MF", "ALL"), run_go_ontology)
names(go_results) <- c("CC", "BP", "MF", "ALL")

message("Done. Outputs written to: ", normalizePath(output_dir, mustWork = FALSE))
