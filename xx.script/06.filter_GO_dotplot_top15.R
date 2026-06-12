############################################################
## 06. GO enrichment stats -> dot-plot top filtering
##
## 목적:
##   03 / 04 / 05 에서 만든 GO_*_enrichment_stats.csv 가 너무 길어서
##   세 데이터셋(GSE37114 / GSE134033 / GSE39872) 비교가 어렵다.
##   각 dot plot 의 cutoff 와 동일한 방식
##     = 지표(metric) × 효과방향(Positive / Negative) 조합마다
##       FDR 가장 유의한 top N term
##   으로 줄여서 비교용 CSV 만 새로 생성한다.
##
## 중요:
##   - plot 을 다시 그리지 않는다. 03/04/05 를 재실행하지 않는다.
##   - 각 폴더 내의 GO_*_enrichment_stats.csv 만 읽어서 작업한다.
##   - dot plot 과 동일 로직: fdr < 0.05 필터 -> direction(effect>=0) ->
##     group_by(metric, direction) %>% slice_min(fdr, n = N).
##     단, plot 은 N = 10 이었고 여기서는 N = 15 로 한다.
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

## 스크립트 위치 기준으로 작업 디렉토리 고정 (03/04/05 와 동일 방식)
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  getwd()
}
setwd(get_script_dir())
cat("Working directory set to:", getwd(), "\n")

base_dir <- "../01.result/own_analysis"

N_TOP      <- 15        # 방향별 상위 개수 (dot plot 은 10 이었음)
ONTOLOGIES <- c("CC", "BP", "MF", "ALL")

## 세 데이터셋. 두 번째 지표 컬럼명이 다르다:
##   GSE37114  : ribosome density change  (ribo_fdr / mean_ribo_density_change)
##   GSE134033 : expression change (KO)   (expr_fdr / mean_expr_change)
##   GSE39872  : expression change (KD)   (expr_fdr / mean_expr_change)
datasets <- list(
  list(
    dir          = "GSE37114_analysis",
    metric2_lab  = "Ribosome density change",
    metric2_fdr  = "ribo_fdr",
    metric2_eff  = "mean_ribo_density_change"
  ),
  list(
    dir          = "GSE134033_analysis",
    metric2_lab  = "Expression change (KO)",
    metric2_fdr  = "expr_fdr",
    metric2_eff  = "mean_expr_change"
  ),
  list(
    dir          = "GSE39872_analysis",
    metric2_lab  = "Expression change (KD)",
    metric2_fdr  = "expr_fdr",
    metric2_eff  = "mean_expr_change"
  )
)

############################################################
## 한 stats CSV -> dot-plot 기준 top N (방향별) long table
############################################################
filter_one <- function(stats_path, ontology, ds) {
  go_stats <- readr::read_csv(stats_path, show_col_types = FALSE)

  ## CLIP enrichment + 두 번째 지표를 long 으로 펼친다 (dot plot 과 동일).
  go_long <- bind_rows(
    go_stats %>% transmute(
      GO, TERM, n_genes,
      metric = "CLIP enrichment",
      fdr    = clip_fdr,
      effect = mean_clip_enrichment
    ),
    go_stats %>% transmute(
      GO, TERM, n_genes,
      metric = ds$metric2_lab,
      fdr    = .data[[ds$metric2_fdr]],
      effect = .data[[ds$metric2_eff]]
    )
  ) %>%
    filter(fdr < 0.05) %>%
    mutate(
      neglog10_fdr = -log10(pmax(fdr, 1e-320)),
      direction    = if_else(effect >= 0, "Positive", "Negative")
    ) %>%
    ## 지표 × 방향 조합마다 FDR 가장 유의한 top N
    group_by(metric, direction) %>%
    slice_min(fdr, n = N_TOP, with_ties = FALSE) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    transmute(
      ontology = ontology,
      metric, direction, rank,
      GO, TERM, n_genes,
      effect, fdr, neglog10_fdr
    ) %>%
    arrange(metric, direction, rank)

  go_long
}

############################################################
## 데이터셋별 실행: 폴더 안에 per-ontology + combined CSV 저장
############################################################
for (ds in datasets) {
  folder <- file.path(base_dir, ds$dir)
  cat("\n== ", ds$dir, " ==\n", sep = "")

  combined <- list()

  for (ont in ONTOLOGIES) {
    stats_path <- file.path(folder, paste0("GO_", ont, "_enrichment_stats.csv"))
    if (!file.exists(stats_path)) {
      cat("  [skip] missing:", basename(stats_path), "\n")
      next
    }

    res <- filter_one(stats_path, ont, ds)

    out_path <- file.path(folder, paste0("GO_", ont, "_dotplot_top", N_TOP, ".csv"))
    readr::write_csv(res, out_path)
    cat("  [ok] ", basename(out_path), " (", nrow(res), " rows)\n", sep = "")

    combined[[ont]] <- res
  }

  ## 네 ontology 를 한 파일로도 합쳐 비교 편하게
  if (length(combined) > 0) {
    all_df <- bind_rows(combined)
    all_path <- file.path(folder, paste0("GO_combined_dotplot_top", N_TOP, ".csv"))
    readr::write_csv(all_df, all_path)
    cat("  [ok] ", basename(all_path), " (", nrow(all_df), " rows)\n", sep = "")
  }
}

cat("\nDone.\n")
