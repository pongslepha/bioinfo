############################################################
## 07. 세 데이터셋 CLIP enrichment GO term 비교 (하나의 dot plot)
##
## 목적:
##   GSE37114 / GSE134033 / GSE39872 세 데이터셋에서 LIN28 "CLIP
##   enrichment" 에 대해 GO term 이 양(Positive) / 음(Negative) 방향으로
##   얼마나 유의하게 enrich 되는지를 하나의 그림으로 비교한다.
##
##   ncc/"03. seurat clustering.R" 의 pathway(ORA) dot plot (라인 540~565)
##   과 동일한 방식:
##     - gotable : 전체 GO 통계표
##     - topn    : 그룹(=여기서는 데이터셋×방향)별 top N term 선택
##     - gotable2: gotable 중 topn 에 뽑힌 term 만, 그리고 유의한 셀만 남김
##     - ggplot  : x=그룹, y=term(Description), size=-log10(p), color=값
##   를 데이터셋 비교용으로 옮긴 것이다.
##
## 입력 (각 폴더의 전체 통계표 csv 만 읽는다. 분석 재실행 없음):
##   ../01.result/own_analysis/<DS>_analysis/GO_ALL_enrichment_stats.csv
##     컬럼: GO, n_genes, mean_clip_enrichment, clip_p, clip_fdr,
##           TERM, ontology, neglog10_fdr ...
##
## 출력 (ontology 별로, 방향 구분 없이 한 그림, GSE37114 순 정렬):
##   ../01.result/own_analysis/comparison_between_datasets/
##     GO_{ALL,CC,BP,MF}_CLIP_enrichment_comparison_top5.png
##     GO_{ALL,CC,BP,MF}_CLIP_enrichment_comparison_top5.csv  (long table)
##
## 실행:
##   ~/miniconda3/envs/lab/bin/Rscript xx.script/07.comparison_between_dataset.R
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

## 스크립트 위치 기준으로 작업 디렉토리 고정 (03/04/05/06 와 동일 방식)
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

base_dir   <- "../01.result/own_analysis"
output_dir <- file.path(base_dir, "comparison_between_datasets")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

N_TOP      <- 5       # 데이터셋 × 방향 별 상위 term 개수
FDR_CUT    <- 0.05
ref_ds     <- "GSE37114"   # y 축 정렬 기준 데이터셋
## ALL(CC+BP+MF 합본) 과 개별 ontology 각각에 대해 그린다
ONTOLOGIES <- c("ALL", "CC", "BP", "MF")

ont_full_name <- c(
  ALL = "All ontologies (CC + BP + MF)",
  CC  = "Cellular Component",
  BP  = "Biological Process",
  MF  = "Molecular Function"
)

## 데이터셋 표시 순서 / 라벨
datasets <- tribble(
  ~dataset,    ~dir,
  "GSE37114",  "GSE37114_analysis",
  "GSE134033", "GSE134033_analysis",
  "GSE39872",  "GSE39872_analysis"
)

############################################################
## 한 ontology 에 대해 gotable -> topn -> gotable2 -> dot plot
############################################################
## seurat 03 의 pathway dot plot 과 동일한 흐름:
##   gotable  : 세 데이터셋의 GO_<ont>_enrichment_stats.csv (CLIP enrichment 부분)
##   topn     : 데이터셋 × 방향 별 top N term (그림은 방향으로 나누지 않음)
##   gotable2 : 선택된 term 을 유의한(fdr<0.05) 모든 데이터셋 셀에 표시
##   plot     : x=데이터셋, y=GO term, size=-log10(FDR), color=CLIP effect(log2)
##              y 축은 GSE37114 의 CLIP enrichment effect 순으로 정렬
run_ontology <- function(ont) {
  cat("\n== GO ", ont, " ==\n", sep = "")

  ## --- 1. gotable: 세 데이터셋 통계표 로드 + CLIP enrichment long ---
  load_one <- function(ds_name, ds_dir) {
    stats_path <- file.path(base_dir, ds_dir,
                            paste0("GO_", ont, "_enrichment_stats.csv"))
    if (!file.exists(stats_path)) {
      warning("missing: ", stats_path)
      return(NULL)
    }
    readr::read_csv(stats_path, show_col_types = FALSE) %>%
      transmute(
        dataset      = ds_name,
        GO, TERM, n_genes,
        effect       = mean_clip_enrichment,
        fdr          = clip_fdr,
        neglog10_fdr = -log10(pmax(clip_fdr, 1e-320)),
        direction    = if_else(mean_clip_enrichment >= 0, "Positive", "Negative")
      )
  }

  gotable <- purrr::pmap_dfr(list(datasets$dataset, datasets$dir), load_one)
  if (nrow(gotable) == 0) {
    message("  [skip] no data for ontology: ", ont)
    return(invisible(NULL))
  }

  ## --- 2. topn: 데이터셋 × 방향 별 top N term (합집합) ---
  topn <- gotable %>%
    filter(fdr < FDR_CUT) %>%
    group_by(dataset, direction) %>%
    slice_min(fdr, n = N_TOP, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(GO, TERM)

  ## --- 3. gotable2: 선택 term 의 유의한 모든 데이터셋 셀 ---
  plot_df <- topn %>%
    left_join(
      gotable %>% dplyr::select(dataset, GO, n_genes, effect, fdr, neglog10_fdr),
      by = "GO", relationship = "many-to-many"
    ) %>%
    filter(fdr < FDR_CUT) %>%
    mutate(dataset = factor(dataset, levels = datasets$dataset))

  cat("  selected terms: ", n_distinct(plot_df$TERM), "\n", sep = "")

  ## y 축(term) 정렬: GSE37114 의 CLIP enrichment effect 순(양성 위, 음성 아래).
  ## GSE37114 에서 유의하지 않은 term(ref_effect=NA)은 맨 아래로.
  term_order <- plot_df %>%
    filter(dataset == ref_ds) %>%
    distinct(TERM, effect) %>%
    dplyr::rename(ref_effect = effect)

  term_levels <- plot_df %>%
    distinct(TERM) %>%
    left_join(term_order, by = "TERM") %>%
    arrange(ref_effect) %>%
    pull(TERM)

  plot_df <- plot_df %>% mutate(TERM = factor(TERM, levels = term_levels))

  ## 그림에 쓰인 long table 저장
  readr::write_csv(
    plot_df %>%
      arrange(dataset, desc(neglog10_fdr)) %>%
      dplyr::select(TERM, GO, dataset, n_genes, effect, fdr, neglog10_fdr),
    file.path(output_dir,
              paste0("GO_", ont, "_CLIP_enrichment_comparison_top", N_TOP, ".csv"))
  )

  ## --- 4. dot plot (seurat 540~565 라인 스타일) ---
  p <- ggplot(plot_df,
              aes(x = dataset, y = TERM,
                  size = neglog10_fdr, color = effect)) +
    geom_point() +
    scale_size_continuous(range = c(1.5, 7), name = "-log10(FDR)") +
    scale_color_gradient2(
      low = "steelblue", mid = "grey90", high = "firebrick",
      midpoint = 0, name = "Mean CLIP\nenrichment (log2)"
    ) +
    labs(
      x = NULL, y = NULL,
      title = paste0("LIN28 CLIP enrichment GO terms across datasets — ",
                     ont_full_name[[ont]]),
      subtitle = paste0("top ", N_TOP, " per direction per dataset; ",
                        "ordered by ", ref_ds,
                        " CLIP enrichment; dot shown where FDR < 0.05")
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major = element_line(color = "grey92"),
      axis.text.x      = element_text(angle = 0, vjust = 1, hjust = 0.5,
                                      face = "bold"),
      axis.text.y      = element_text(size = 7),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5, size = 7),
      legend.position  = "right",
      legend.title     = element_text(size = 8)
    )

  n_terms <- n_distinct(plot_df$TERM)
  plot_h  <- max(5, 0.2 * n_terms + 2)

  ggsave(
    file.path(output_dir,
              paste0("GO_", ont, "_CLIP_enrichment_comparison_top", N_TOP, ".png")),
    p, width = 9, height = plot_h, dpi = 900, units = "in", limitsize = FALSE
  )
}

for (ont in ONTOLOGIES) run_ontology(ont)

cat("\nDone. Outputs written to:",
    normalizePath(output_dir, mustWork = FALSE), "\n")
