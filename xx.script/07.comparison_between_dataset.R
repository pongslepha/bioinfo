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
## 출력 (Positive / Negative 를 한 그림에 위/아래 패널로):
##   ../01.result/own_analysis/comparison_between_datasets/
##     GO_CLIP_enrichment_comparison_top5.png
##     GO_CLIP_enrichment_comparison_top5.csv  (그림에 쓰인 long table)
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

N_TOP    <- 5         # 데이터셋 × 방향 별 상위 term 개수
ONTOLOGY <- "ALL"     # CC+BP+MF 합친 통계표 사용 (03/04/05 의 ALL 과 동일)
FDR_CUT  <- 0.05

## 데이터셋 표시 순서 / 라벨
datasets <- tribble(
  ~dataset,    ~dir,
  "GSE37114",  "GSE37114_analysis",
  "GSE134033", "GSE134033_analysis",
  "GSE39872",  "GSE39872_analysis"
)

############################################################
## 1. 각 데이터셋의 전체 GO 통계표(gotable) 로드 + CLIP enrichment 부분만
############################################################
## seurat 의 gotable 에 해당. 여기서는 "CLIP enrichment" 지표만 비교하므로
## clip_fdr / mean_clip_enrichment 만 long 형태로 가져온다.
load_one <- function(ds_name, ds_dir) {
  stats_path <- file.path(base_dir, ds_dir,
                          paste0("GO_", ONTOLOGY, "_enrichment_stats.csv"))
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

gotable <- purrr::pmap_dfr(
  list(datasets$dataset, datasets$dir),
  load_one
)

stopifnot(nrow(gotable) > 0)
cat("Loaded GO stats rows:", nrow(gotable),
    "from", n_distinct(gotable$dataset), "datasets\n")

############################################################
## 2. topn : 데이터셋 × 방향 별 top N term 선택 (seurat 의 topn)
############################################################
## 각 데이터셋에서 CLIP enrichment 가 가장 유의한 term 을 Positive / Negative
## 방향마다 N_TOP 개씩 뽑는다. 이 term 들이 그림의 y 축 후보가 된다.
topn <- gotable %>%
  filter(fdr < FDR_CUT) %>%
  group_by(dataset, direction) %>%
  slice_min(fdr, n = N_TOP, with_ties = FALSE) %>%
  ungroup() %>%
  ## term 은 어느 한 방향으로 선택되면 그 방향 패널에 들어간다.
  distinct(direction, GO, TERM)

cat("Selected (direction, GO) pairs:", nrow(topn), "\n")

############################################################
## 3. gotable2 : 선택된 term 을, 유의한(=fdr<0.05) 모든 데이터셋 셀에서
##              표시 (seurat 의 gotable2)
############################################################
## 한 term 이 GSE39872 에서만 top5 이었더라도, 다른 데이터셋에서 유의하면
## 그 데이터셋 칸에도 점을 찍어 데이터셋 간 비교가 되게 한다.
plot_df <- topn %>%
  ## facet 은 "선택된 방향" 기준 (Positive / Negative 패널)
  dplyr::rename(facet_direction = direction) %>%
  left_join(
    gotable %>% dplyr::select(dataset, GO, n_genes, effect, fdr, neglog10_fdr),
    by = "GO", relationship = "many-to-many"
  ) %>%
  filter(fdr < FDR_CUT) %>%
  mutate(
    dataset         = factor(dataset, levels = datasets$dataset),
    facet_direction = factor(facet_direction, levels = c("Positive", "Negative"))
  )

## y 축(term) 정렬: 패널마다 독립적으로 가장 강한 신호(neglog10_fdr) 순.
## 같은 term 이 Positive / Negative 두 패널에 모두 나올 수 있으므로
## "facet||TERM" 을 factor key 로 써서 패널별 정렬을 보장한다.
term_order <- plot_df %>%
  group_by(facet_direction, TERM) %>%
  summarise(ord_val = max(neglog10_fdr), .groups = "drop") %>%
  arrange(facet_direction, ord_val)

plot_df <- plot_df %>%
  mutate(
    y_key = factor(
      paste(facet_direction, TERM, sep = "||"),
      levels = paste(term_order$facet_direction, term_order$TERM, sep = "||")
    )
  )

## 그림에 쓰인 long table 저장
readr::write_csv(
  plot_df %>%
    arrange(facet_direction, desc(neglog10_fdr), dataset) %>%
    dplyr::select(facet_direction, TERM, GO, dataset,
                  n_genes, effect, fdr, neglog10_fdr),
  file.path(output_dir, "GO_CLIP_enrichment_comparison_top5.csv")
)

############################################################
## 4. 통합 dot plot (seurat 540~565 라인 스타일)
##    x = 데이터셋, y = GO term, size = -log10(FDR), color = CLIP effect(log2)
############################################################
p <- ggplot(plot_df,
            aes(x = dataset, y = y_key,
                size = neglog10_fdr, color = effect)) +
  geom_point() +
  facet_grid(facet_direction ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
  scale_size_continuous(range = c(1.5, 7), name = "-log10(FDR)") +
  ## 효과 방향이 보이도록 발산형(diverging) 팔레트, 0 을 흰색으로
  scale_color_gradient2(
    low = "steelblue", mid = "grey90", high = "firebrick",
    midpoint = 0, name = "Mean CLIP\nenrichment (log2)"
  ) +
  labs(
    x = NULL, y = NULL,
    title = paste0("LIN28 CLIP enrichment GO terms across datasets",
                   " (top ", N_TOP, " per direction per dataset)"),
    subtitle = "GO_ALL_enrichment_stats; dot shown where FDR < 0.05"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major   = element_line(color = "grey92"),
    axis.text.x        = element_text(angle = 0, vjust = 1, hjust = 0.5,
                                       face = "bold"),
    axis.text.y        = element_text(size = 7),
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(hjust = 0.5, size = 7),
    legend.position    = "right",
    legend.title       = element_text(size = 8)
  )

## 패널별 term 수에 맞춰 세로 길이 자동 조정
n_terms <- n_distinct(plot_df$y_key)
plot_h  <- max(5, 0.2 * n_terms + 2)

ggsave(
  file.path(output_dir, "GO_CLIP_enrichment_comparison_top5.png"),
  p, width = 9, height = plot_h, dpi = 900, units = "in", limitsize = FALSE
)

cat("Done. Outputs written to:",
    normalizePath(output_dir, mustWork = FALSE), "\n")
