############################################################
## GSE134033 LIN28A CLIP + RNA-seq (knockout) pathway analysis
##
## GSE134033 ("LIN28A binds to meiotic gene transcripts and
## modulates their translation in male germ cells", mouse testes)
## 은 원본 03.pathway_analysis.R 과 동일한 LIN28A 를 다룬다.
## 따라서 동일한 틀(CLIP enrichment vs. 기능적 변화 + GO Wilcoxon
## enrichment + bubble plot)을 그대로 적용한다.
##
## 원본과 다른 점:
##   - RPF(ribosome profiling)가 없으므로 y축을 "ribosome density
##     change" 대신 RNA-seq 의 "Lin28a KO vs control 발현 변화"로 둔다.
##   - CLIP / RNA 가 별도 파일이라 두 metric 을 직접 계산해 join 한다.
##
## Input:
##   00.data/GEO/GSE134033_CLIP-seq_gene_expression.xlsx
##       CLIP-seq(LIN28A IP) FPKM: DP5i1/DP5i3/DP5i7 (triplicate)
##       gene_id = ENSMUSG, gene_short_name = symbol
##   00.data/GEO/GSE134033_10d_RNA-seq_genes_fpkm_annotation.csv.gz
##       cuffdiff 결과(10dpp testes): case = Lin28a KO, control = WT
##       gene_id = RefSeq(NM_*), geneName = symbol
##   00.data/mouselocalization-20210507.txt   (optional, ENSMUSG 로 join)
############################################################

library(tidyverse)
library(data.table)
library(readxl)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(ggrepel)

############################################################
## 0. Paths  (xx.script/ 에서 실행한다고 가정 - 03 스크립트와 동일)
############################################################

data_dir   <- "../00.data"
geo_dir    <- file.path(data_dir, "GEO")
output_dir <- "../01.result/own_analysis/GSE134033_analysis"

clip_file <- file.path(geo_dir, "GSE134033_CLIP-seq_gene_expression.xlsx")
rna_file  <- file.path(geo_dir, "GSE134033_10d_RNA-seq_genes_fpkm_annotation.csv.gz")
loc_file  <- file.path(data_dir, "mouselocalization-20210507.txt")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## case = Lin28a KO, control = WT 로 가정한다 (GEO design 기준).
## y축 부호 해석이 반대라면 이 한 줄의 라벨/계산만 뒤집으면 된다.
pc <- 1  # log 계산용 pseudocount (원본과 동일)

############################################################
## 1. Load CLIP-seq (LIN28A IP) FPKM
############################################################
## DP5i1/DP5i3/DP5i7 세 replicate 의 FPKM 평균을 CLIP 신호로 사용.

clip_raw <- read_excel(clip_file, sheet = "gene_expression")

clip_df <- clip_raw %>%
  transmute(
    ensembl_gene_id = sub("\\..*$", "", gene_id),
    symbol          = toupper(gene_short_name),
    clip_fpkm_mean  = (DP5i1_FPKM + DP5i3_FPKM + DP5i7_FPKM) / 3
  ) %>%
  filter(ensembl_gene_id != "", !is.na(symbol)) %>%
  ## ENSMUSG 당 한 행 (혹시 모를 중복 대비, 신호 큰 쪽 유지)
  group_by(ensembl_gene_id) %>%
  slice_max(clip_fpkm_mean, n = 1, with_ties = FALSE) %>%
  ungroup()

############################################################
## 2. Load RNA-seq (cuffdiff: Lin28a KO vs control) FPKM
############################################################
## RefSeq 기반이라 symbol 로 join 한다. symbol 당 가장 발현 높은
## isoform 한 행만 남긴다.

rna_raw <- fread(cmd = paste("zcat", shQuote(rna_file)))

rna_df <- rna_raw %>%
  as_tibble() %>%
  transmute(
    symbol      = toupper(geneName),
    rna_case    = suppressWarnings(as.numeric(case)),     # Lin28a KO FPKM
    rna_control = suppressWarnings(as.numeric(control))   # WT FPKM
  ) %>%
  filter(symbol != "-", symbol != "", !is.na(symbol),
         !is.na(rna_case), !is.na(rna_control)) %>%
  group_by(symbol) %>%
  slice_max(rna_control, n = 1, with_ties = FALSE) %>%
  ungroup()

############################################################
## 3. Join CLIP + RNA and compute the two metrics
############################################################
## clip_enrichment : LIN28A 가 얼마나 그 전사체에 결합하는가
##                   = log2(CLIP IP FPKM / WT RNA FPKM)   (원본 x축과 동일 개념)
## expr_change     : Lin28a KO 시 mRNA 수준 변화
##                   = log2(KO FPKM / WT FPKM)            (원본 y축 자리 대체)

gene_stat <- clip_df %>%
  inner_join(rna_df, by = "symbol") %>%
  mutate(
    clip_enrichment = log2((clip_fpkm_mean + pc) / (rna_control + pc)),
    expr_change     = log2((rna_case + pc) / (rna_control + pc))
  )

############################################################
## 4. Basic filtering
############################################################
## 발현이 거의 없는 유전자는 제거하고, 원본과 동일하게 histone 유전자
## (poly(A) 정규화 artifact 가능성) 도 제거한다.

gene_stat <- gene_stat %>%
  filter(
    (rna_control + rna_case) >= 1,
    is.finite(clip_enrichment),
    is.finite(expr_change),
    !str_detect(symbol, "^HIST")
  )

message(sprintf("Genes after filtering: %d", nrow(gene_stat)))

############################################################
## 5. Load mouse localization file and join (ENSMUSG)
############################################################
## localization 파일은 ENSMUSG 기준이라 CLIP gene_id 와 직접 join 가능.

loc_raw <- fread(loc_file)

loc_df <- loc_raw %>%
  as_tibble() %>%
  transmute(
    ensembl_gene_id    = sub("\\..*$", "", as.character(gene_id)),
    localization_group = str_to_title(as.character(type))
  ) %>%
  filter(ensembl_gene_id != "", localization_group != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

gene_stat_loc <- gene_stat %>%
  left_join(loc_df, by = "ensembl_gene_id") %>%
  mutate(localization_group = replace_na(localization_group, "Unannotated"))

############################################################
## 6. Save gene-level table
############################################################

write_csv(
  gene_stat_loc,
  file.path(output_dir, "gene_level_CLIP_RNA_localization_table.csv")
)

############################################################
## 7. Localization group-level summary + Wilcoxon test
############################################################
## 원본 11~13 단계와 동일: 각 localization group vs 나머지.

loc_summary <- gene_stat_loc %>%
  filter(localization_group != "Unannotated") %>%
  group_by(localization_group) %>%
  summarise(
    n_genes                  = n_distinct(ensembl_gene_id),
    mean_clip_enrichment     = mean(clip_enrichment, na.rm = TRUE),
    median_clip_enrichment   = median(clip_enrichment, na.rm = TRUE),
    mean_expr_change         = mean(expr_change, na.rm = TRUE),
    median_expr_change       = median(expr_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_clip_enrichment))

all_loc_genes <- gene_stat_loc %>%
  filter(localization_group != "Unannotated") %>%
  dplyr::select(ensembl_gene_id, symbol, localization_group,
                clip_enrichment, expr_change) %>%
  distinct()

test_one_loc <- function(loc_group) {
  in_df  <- all_loc_genes %>% filter(localization_group == loc_group)
  out_df <- all_loc_genes %>% filter(localization_group != loc_group)
  if (nrow(in_df) < 10 || nrow(out_df) < 10) return(NULL)

  tibble(
    localization_group   = loc_group,
    n_genes              = nrow(in_df),
    mean_clip_enrichment = mean(in_df$clip_enrichment, na.rm = TRUE),
    mean_expr_change     = mean(in_df$expr_change, na.rm = TRUE),
    clip_p = wilcox.test(in_df$clip_enrichment, out_df$clip_enrichment)$p.value,
    expr_p = wilcox.test(in_df$expr_change,     out_df$expr_change)$p.value
  )
}

loc_test <- map_dfr(unique(all_loc_genes$localization_group), test_one_loc) %>%
  mutate(
    clip_fdr     = p.adjust(clip_p, method = "BH"),
    expr_fdr     = p.adjust(expr_p, method = "BH"),
    combined_fdr = pmin(clip_fdr, expr_fdr, na.rm = TRUE),
    neglog10_fdr = -log10(combined_fdr)
  ) %>%
  arrange(combined_fdr)

loc_group_stats <- loc_test %>%
  left_join(
    loc_summary %>%
      dplyr::select(localization_group, median_clip_enrichment, median_expr_change),
    by = "localization_group"
  )

write_csv(
  loc_group_stats,
  file.path(output_dir, "localization_group_summary_and_test.csv")
)

print(loc_group_stats)

############################################################
## 8. Localization group bubble plot (원본 13 단계와 동일)
############################################################

p_loc_bubble <- ggplot(
  loc_test,
  aes(x = mean_clip_enrichment, y = mean_expr_change)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(size = n_genes, color = neglog10_fdr), alpha = 0.8) +
  geom_text_repel(
    aes(label = paste0(localization_group, " (", n_genes, ")")),
    size = 4, max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(4, 18), name = "Number of genes") +
  scale_color_gradient(low = "khaki1", high = "red4", name = "-log10(FDR)") +
  labs(
    title = "GSE134033: localization-group enrichment (CLIP vs Lin28a-KO expression)",
    x = "Mean LIN28A CLIP enrichment (log2)",
    y = "Mean log2 fold change (Lin28a KO / control)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_line(color = "grey85", linetype = "dashed"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "localization_group_bubble_plot.png"),
  p_loc_bubble, width = 9, height = 6, dpi = 900
)

############################################################
## 9. GO enrichment Figure 5A-like plots (CC / BP / MF)
############################################################
## 원본 14 단계와 동일. ENSMUSG -> ENTREZ/SYMBOL -> GO 매핑 후
## ontology 별 group-wise Wilcoxon -> 통계표 + 2D bubble plot.

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

go_map_all <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(gene_stat_go$ENTREZID),
  keytype = "ENTREZID",
  columns = c("GO", "ONTOLOGY")
) %>%
  filter(!is.na(GO)) %>%
  distinct(ENTREZID, GO, .keep_all = TRUE)

all_genes_go <- gene_stat_go %>%
  dplyr::select(ENTREZID, clip_enrichment, expr_change) %>%
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

  dat_go <- go_map %>% inner_join(gene_stat_go, by = "ENTREZID")

  test_one_go <- function(go_id) {
    genes_in <- dat_go %>% filter(GO == go_id) %>% pull(ENTREZID) %>% unique()
    in_df  <- all_genes_go %>% filter(ENTREZID %in% genes_in)
    out_df <- all_genes_go %>% filter(!ENTREZID %in% genes_in)
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

  label_go <- go_stats %>% slice_head(n = 25)

  p_go <- ggplot(
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
      title = paste0("GSE134033 GO ", ont_full_name[[ont]],
                     " analysis (LIN28A CLIP vs Lin28a-KO expression)"),
      x = "Mean LIN28A CLIP enrichment (log2)",
      y = "Mean log2 fold change (Lin28a KO / control)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_line(color = "grey85", linetype = "dashed"),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(output_dir, paste0("GO_", ont, "_Figure5A_like.png")),
    p_go, width = 12, height = 7, dpi = 900
  )

  ## ----- Dot (lollipop) plot ----- (원본 03 스크립트와 동일한 스타일)
  ## 두 지표(CLIP enrichment / Lin28a-KO expression change) × 효과 방향
  ## (Positive/Negative) 으로 facet 을 나누고, 각 칸마다 FDR 상위 n_top 개
  ## term 을 -log10(FDR) lollipop 으로 표시. 점 크기 = 유전자 수.
  n_top <- 10

  ## facet_wrap(ncol = 2) 은 panel 을 행 우선(row-major)으로 채우므로,
  ## level 순서를 row1=Negative, row2=Positive / col1=CLIP, col2=Expression 으로 맞춘다.
  metric_levels <- c("CLIP enrichment", "Expression change (KO)")  # column 순서
  dir_levels    <- c("Negative", "Positive")                       # row 순서

  go_long <- bind_rows(
    go_stats %>% transmute(
      TERM, n_genes, metric = "CLIP enrichment",
      fdr = clip_fdr, effect = mean_clip_enrichment
    ),
    go_stats %>% transmute(
      TERM, n_genes, metric = "Expression change (KO)",
      fdr = expr_fdr, effect = mean_expr_change
    )
  ) %>%
    filter(fdr < 0.05) %>%
    mutate(
      neglog10_fdr = -log10(pmax(fdr, 1e-320)),  # FDR==0 -> 320 으로 cap
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
        "GSE134033 GO ", ont_full_name[[ont]],
        " enrichment (LIN28A CLIP vs Lin28a-KO expression; top ", n_top,
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
