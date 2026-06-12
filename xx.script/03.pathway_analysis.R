############################################################
## LIN28A CLIP + RPF + localization analysis
## Input:
##   01.result/read-counts.txt
##   00.data/mouselocalization-20210507.txt
##   00.data/gencode.vM27.annotation.gtf.gz
############################################################

library(tidyverse)
library(data.table)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(ggrepel)

############################################################
## 0. Paths
############################################################
## 상대경로(../01.result 등)는 작업 디렉토리 기준이라, 실행 위치에
## 따라 깨진다. 스크립트(.R)가 있는 폴더(xx.script)를 찾아 그 기준으로
## 작업 디렉토리를 맞춰, 어디서 실행하든 동일하게 동작하게 한다.
get_script_dir <- function() {
  ## 1) Rscript 실행: --file= 인자
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  ## 2) source("...") 실행: sys.frame 의 ofile
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  ## 3) RStudio "Source" 버튼 등: 활성 문서 경로
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  ## 4) 모두 실패하면 현재 작업 디렉토리 사용
  getwd()
}

setwd(get_script_dir())
cat("Working directory set to:", getwd(), "\n")

result_dir <- "../01.result"
data_dir   <- "../00.data"
output_dir <- "../01.result/own_analysis/GSE37114_analysis"

read_count_file <- file.path(result_dir, "read-counts.txt")
gtf_file        <- file.path(data_dir, "gencode.vM27.annotation.gtf.gz")
loc_file        <- file.path(data_dir, "mouselocalization-20210507.txt")

## 결과물(csv/plot)은 output_dir 에 저장
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Load featureCounts result
############################################################

fc <- fread(read_count_file, skip = 1)

## featureCounts는 BAM 컬럼명을 전체 경로
## (예: "../00.data/extracted/binfo1-work/CLIP-35L33G.bam") 로 기록하므로
## basename 으로 정리해서 아래 select 에서 짧은 이름으로 참조 가능하게 함.
bam_cols <- grepl("\\.bam$", colnames(fc))
colnames(fc)[bam_cols] <- basename(colnames(fc)[bam_cols])

fc <- fc %>%
  mutate(
    ensembl_gene_id_version = Geneid,
    ensembl_gene_id = sub("\\..*$", "", Geneid)
  )

count_df <- fc %>%
  dplyr::select(
    ensembl_gene_id,
    Length,
    CLIP = `CLIP-35L33G.bam`,
    RNA_control = `RNA-control.bam`,
    RNA_siLin28a = `RNA-siLin28a.bam`,
    RNA_siLuc = `RNA-siLuc.bam`,
    RPF_siLin28a = `RPF-siLin28a.bam`,
    RPF_siLuc = `RPF-siLuc.bam`
  )

############################################################
## 2. Basic filtering
############################################################

count_df <- count_df %>%
  mutate(
    total_RNA  = RNA_control + RNA_siLin28a + RNA_siLuc,
    total_RPF  = RPF_siLin28a + RPF_siLuc,
    total_CLIP = CLIP
  ) %>%
  filter(
    total_RNA >= 20,
    total_RPF >= 10
  )

############################################################
## 3. Library-size normalization and metric calculation
############################################################

cpm <- function(x) {
  x / sum(x, na.rm = TRUE) * 1e6
}

pc <- 1

count_df <- count_df %>%
  mutate(
    CLIP_cpm = cpm(CLIP),
    RNA_control_cpm = cpm(RNA_control),
    RNA_siLin28a_cpm = cpm(RNA_siLin28a),
    RNA_siLuc_cpm = cpm(RNA_siLuc),
    RPF_siLin28a_cpm = cpm(RPF_siLin28a),
    RPF_siLuc_cpm = cpm(RPF_siLuc)
  ) %>%
  mutate(
    clip_enrichment =
      log2((CLIP_cpm + pc) / (RNA_control_cpm + pc)),

    ribo_density_siLin28a =
      log2((RPF_siLin28a_cpm + pc) / (RNA_siLin28a_cpm + pc)),

    ribo_density_siLuc =
      log2((RPF_siLuc_cpm + pc) / (RNA_siLuc_cpm + pc)),

    ribo_density_change =
      ribo_density_siLin28a - ribo_density_siLuc
  )

############################################################
## 4. Add gene symbol from GENCODE GTF
############################################################
## mouselocalization file이 gene symbol 기준일 가능성이 높으므로
## GTF에서 gene_id - gene_name mapping을 가져옴.
############################################################

## rtracklayer 없이 GTF 의 gene line 만 직접 파싱 (gene_id / gene_name / gene_type).
gtf_raw <- fread(
  cmd = paste("zcat", shQuote(gtf_file)),
  sep = "\t", header = FALSE, quote = "",
  col.names = c("seqname", "source", "feature", "start", "end",
                "score", "strand", "frame", "attribute")
)

gtf_gene <- gtf_raw %>%
  filter(feature == "gene") %>%
  transmute(
    ensembl_gene_id = sub("\\..*$", "",
                          sub('.*gene_id "([^"]+)".*', "\\1", attribute)),
    gene_name = sub('.*gene_name "([^"]+)".*', "\\1", attribute),
    gene_type = sub('.*gene_type "([^"]+)".*', "\\1", attribute)
  ) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

gene_stat <- count_df %>%
  left_join(gtf_gene, by = "ensembl_gene_id")

############################################################
## 5. Remove problematic genes if needed
############################################################
## Histone genes는 poly(A) RNA-seq normalization artifact 가능성이 있음.
############################################################

gene_stat <- gene_stat %>%
  filter(
    !is.na(gene_name),
    !str_detect(gene_name, "^Hist")
  )

############################################################
## 6. Load mouse localization file
############################################################

loc_raw <- fread(loc_file)

############################################################
## 7. Build localization table (gene_id / type)
############################################################
## 파일 구조: gene_id(ENSEMBL) / "Gene names" / type(=localization).
## ENSEMBL ID 로 직접 join 한다 (symbol 매칭보다 정확).
############################################################

loc_df <- loc_raw %>%
  as_tibble() %>%
  transmute(
    ensembl_gene_id = sub("\\..*$", "", as.character(gene_id)),
    localization_group = str_to_title(as.character(type))
  ) %>%
  filter(ensembl_gene_id != "", localization_group != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

############################################################
## 9. Join localization to gene-level statistics
############################################################

gene_stat_loc <- gene_stat %>%
  left_join(loc_df, by = "ensembl_gene_id") %>%
  mutate(
    localization_group = replace_na(localization_group, "Unannotated")
  )

############################################################
## 10. Save gene-level table
############################################################

write_csv(
  gene_stat_loc,
  file.path(output_dir, "gene_level_CLIP_RPF_localization_table.csv")
)

############################################################
## 11. Localization group-level summary
############################################################

loc_summary <- gene_stat_loc %>%
  filter(localization_group != "Unannotated") %>%
  group_by(localization_group) %>%
  summarise(
    n_genes = n_distinct(ensembl_gene_id),
    mean_clip_enrichment = mean(clip_enrichment, na.rm = TRUE),
    median_clip_enrichment = median(clip_enrichment, na.rm = TRUE),
    mean_ribo_density_change = mean(ribo_density_change, na.rm = TRUE),
    median_ribo_density_change = median(ribo_density_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_clip_enrichment))

############################################################
## 12. Statistical test for each localization group
############################################################
## 각 localization group 유전자 vs 나머지 유전자
############################################################

all_loc_genes <- gene_stat_loc %>%
  filter(localization_group != "Unannotated") %>%
  dplyr::select(
    ensembl_gene_id,
    gene_name,
    localization_group,
    clip_enrichment,
    ribo_density_change
  ) %>%
  distinct()

test_one_loc <- function(loc_group) {
  in_df <- all_loc_genes %>%
    filter(localization_group == loc_group)

  out_df <- all_loc_genes %>%
    filter(localization_group != loc_group)

  if (nrow(in_df) < 10 || nrow(out_df) < 10) {
    return(NULL)
  }

  tibble(
    localization_group = loc_group,
    n_genes = nrow(in_df),
    mean_clip_enrichment = mean(in_df$clip_enrichment, na.rm = TRUE),
    mean_ribo_density_change = mean(in_df$ribo_density_change, na.rm = TRUE),
    clip_p = wilcox.test(
      in_df$clip_enrichment,
      out_df$clip_enrichment
    )$p.value,
    ribo_p = wilcox.test(
      in_df$ribo_density_change,
      out_df$ribo_density_change
    )$p.value
  )
}

loc_test <- map_dfr(
  unique(all_loc_genes$localization_group),
  test_one_loc
) %>%
  mutate(
    clip_fdr = p.adjust(clip_p, method = "BH"),
    ribo_fdr = p.adjust(ribo_p, method = "BH"),
    combined_fdr = pmin(clip_fdr, ribo_fdr, na.rm = TRUE),
    neglog10_fdr = -log10(combined_fdr)
  ) %>%
  arrange(combined_fdr)

## summary(평균/중앙값)와 통계검정(p/FDR) 결과를 한 파일로 합쳐 저장
loc_group_stats <- loc_test %>%
  left_join(
    loc_summary %>%
      dplyr::select(
        localization_group,
        median_clip_enrichment,
        median_ribo_density_change
      ),
    by = "localization_group"
  )

write_csv(
  loc_group_stats,
  file.path(output_dir, "localization_group_summary_and_test.csv")
)

print(loc_group_stats)

############################################################
## 13. Localization group bubble plot
############################################################
## localization group 자체를 GO term처럼 하나의 bubble로 표시
############################################################

p_loc_bubble <- ggplot(
  loc_test,
  aes(
    x = mean_clip_enrichment,
    y = mean_ribo_density_change
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(
    aes(
      size = n_genes,
      color = neglog10_fdr
    ),
    alpha = 0.8
  ) +
  geom_text_repel(
    aes(label = paste0(localization_group, " (", n_genes, ")")),
    size = 4,
    max.overlaps = Inf
  ) +
  scale_size_continuous(
    range = c(4, 18),
    name = "Number of genes"
  ) +
  scale_color_gradient(
    low = "khaki1",
    high = "red4",
    name = "-log10(FDR)"
  ) +
  labs(
    title = "Localization-group enrichment for CLIP and ribosome profiling",
    x = "Mean LIN28A CLIP enrichment (log2)",
    y = "Mean ribosome density change upon Lin28a knockdown (log2)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_line(color = "grey85", linetype = "dashed"),
    plot.title = element_text(face = "bold")
  )

print(p_loc_bubble)

ggsave(
  file.path(output_dir, "localization_group_bubble_plot.png"),
  p_loc_bubble,
  width = 9,
  height = 6,
  dpi = 900
)

############################################################
## 14. GO enrichment Figure 5A-like plots (CC / BP / MF)
############################################################
## localization 파일과 별개로 GO Cellular Component / Biological
## Process / Molecular Function 각 ontology 별로 동일한 분석/plot 생성.
############################################################

gene_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = gene_stat$ensembl_gene_id,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

gene_stat_go <- gene_stat %>%
  left_join(gene_map, by = c("ensembl_gene_id" = "ENSEMBL")) %>%
  filter(!is.na(ENTREZID))

## GO 매핑은 세 ontology 공통이므로 한 번만 가져온 뒤 ontology 별로 나눈다.
go_map_all <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(gene_stat_go$ENTREZID),
  keytype = "ENTREZID",
  columns = c("GO", "ONTOLOGY")
) %>%
  filter(!is.na(GO)) %>%
  distinct(ENTREZID, GO, .keep_all = TRUE)

all_genes_go <- gene_stat_go %>%
  dplyr::select(ENTREZID, clip_enrichment, ribo_density_change) %>%
  distinct()

ont_full_name <- c(
  CC = "Cellular Component",
  BP = "Biological Process",
  MF = "Molecular Function",
  ALL = "All ontologies (CC + BP + MF)"
)

## 한 ontology 에 대해 group-wise Wilcoxon 검정 -> 통계표 + bubble plot 저장
## ont == "ALL" 이면 세 ontology 를 구분 없이 합쳐서 분석한다.
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

  dat_go <- go_map %>%
    inner_join(gene_stat_go, by = "ENTREZID")

  test_one_go <- function(go_id) {
    genes_in <- dat_go %>%
      filter(GO == go_id) %>%
      pull(ENTREZID) %>%
      unique()

    in_df <- all_genes_go %>% filter(ENTREZID %in% genes_in)
    out_df <- all_genes_go %>% filter(!ENTREZID %in% genes_in)

    if (nrow(in_df) < 10 || nrow(out_df) < 10) {
      return(NULL)
    }

    tibble(
      GO = go_id,
      n_genes = nrow(in_df),
      mean_clip_enrichment = mean(in_df$clip_enrichment, na.rm = TRUE),
      mean_ribo_density_change = mean(in_df$ribo_density_change, na.rm = TRUE),
      clip_p = wilcox.test(
        in_df$clip_enrichment,
        out_df$clip_enrichment
      )$p.value,
      ribo_p = wilcox.test(
        in_df$ribo_density_change,
        out_df$ribo_density_change
      )$p.value
    )
  }

  go_stats <- map_dfr(unique(dat_go$GO), test_one_go) %>%
    left_join(go_term_names, by = c("GO" = "GOID")) %>%
    mutate(
      ontology = ont,
      clip_fdr = p.adjust(clip_p, method = "BH"),
      ribo_fdr = p.adjust(ribo_p, method = "BH"),
      combined_fdr = pmin(clip_fdr, ribo_fdr, na.rm = TRUE),
      neglog10_fdr = -log10(combined_fdr)
    ) %>%
    filter(
      n_genes >= 10,
      n_genes <= 500,
      combined_fdr < 0.05
    ) %>%
    arrange(combined_fdr)

  write_csv(
    go_stats,
    file.path(output_dir, paste0("GO_", ont, "_enrichment_stats.csv"))
  )

  ## ----- Figure 5A-like bubble plot -----
  ## term 마다 점 하나: x = mean CLIP enrichment, y = mean ribosome density
  ## change, 점 크기 = 유전자 수, 색 = -log10(FDR). 유의도 상위 25개만 라벨링.
  label_go <- go_stats %>% slice_head(n = 25)

  p_bubble <- ggplot(
    go_stats %>% arrange(neglog10_fdr),  # 유의한(빨간) 점이 맨 위에 그려지도록
    aes(
      x = mean_clip_enrichment,
      y = mean_ribo_density_change
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(
      aes(size = n_genes, color = neglog10_fdr),
      alpha = 0.75
    ) +
    geom_text_repel(
      data = label_go,
      aes(label = paste0(TERM, " (", n_genes, ")")),
      size = 3,
      max.overlaps = Inf
    ) +
    scale_size_continuous(range = c(2, 14), name = "Number of genes") +
    scale_color_gradient(
      low = "khaki1", high = "red4", name = "-log10(FDR)"
    ) +
    labs(
      title = paste0(
        "GO ", ont_full_name[[ont]],
        " analysis for CLIP and ribosome profiling"
      ),
      x = "Mean LIN28A CLIP enrichment (log2)",
      y = "Mean ribosome density change upon Lin28a knockdown (log2)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_line(color = "grey85", linetype = "dashed"),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(output_dir, paste0("GO_", ont, "_Figure5A_like.png")),
    p_bubble,
    width = 12,
    height = 7,
    dpi = 900
  )

  ## ----- Dot (lollipop) plot -----
  ## 두 지표(CLIP enrichment / ribosome density change) × 효과 방향(Positive/
  ## Negative) 으로 facet 을 나누고, 각 칸마다 FDR 상위 n_top 개 term 을
  ## -log10(FDR) lollipop 으로 표시. 점 크기 = 유전자 수.
  n_top <- 10

  ## facet_wrap(ncol = 2) 은 panel 을 행 우선(row-major)으로 채우므로,
  ## level 순서를 row1=Negative, row2=Positive / col1=CLIP, col2=Ribosome 으로 맞춘다.
  metric_levels <- c("CLIP enrichment", "Ribosome density change")  # column 순서
  dir_levels    <- c("Negative", "Positive")                        # row 순서

  go_long <- bind_rows(
    go_stats %>% transmute(
      TERM, n_genes, metric = "CLIP enrichment",
      fdr = clip_fdr, effect = mean_clip_enrichment
    ),
    go_stats %>% transmute(
      TERM, n_genes, metric = "Ribosome density change",
      fdr = ribo_fdr, effect = mean_ribo_density_change
    )
  ) %>%
    filter(fdr < 0.05) %>%
    mutate(
      neglog10_fdr = -log10(pmax(fdr, 1e-320)),  # FDR==0 -> 320 으로 cap
      direction = if_else(effect >= 0, "Positive", "Negative")
    ) %>%
    ## 지표 × 방향 조합마다 가장 유의한 top n_top term
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
    ## facet 별로 term 정렬이 독립적이도록 row 기반 factor 사용
    arrange(facet_lab, neglog10_fdr) %>%
    mutate(term_ord = factor(row_number(), levels = row_number(), labels = TERM))

  p_go <- ggplot(
    go_long,
    aes(x = neglog10_fdr, y = term_ord)
  ) +
    geom_segment(
      aes(x = 0, xend = neglog10_fdr, yend = term_ord, color = direction)
    ) +
    geom_point(aes(size = n_genes, color = direction)) +
    facet_wrap(~ facet_lab, scales = "free_y", ncol = 2) +
    scale_color_manual(
      values = c(Positive = "firebrick", Negative = "steelblue"),
      name = "Effect direction"
    ) +
    scale_size_continuous(range = c(2, 9), name = "Number of genes") +
    labs(
      title = paste0(
        "GO ", ont_full_name[[ont]],
        " enrichment for CLIP and ribosome profiling (top ", n_top,
        " per direction)"
      ),
      x = "-log10(FDR)",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(output_dir, paste0("GO_", ont, "_dotplot.png")),
    p_go,
    width = 16,
    height = 8,
    dpi = 900
  )

  go_stats
}

go_results <- lapply(c("CC", "BP", "MF", "ALL"), run_go_ontology)
names(go_results) <- c("CC", "BP", "MF", "ALL")