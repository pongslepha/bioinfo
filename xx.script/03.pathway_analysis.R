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

result_dir <- "../01.result"
data_dir   <- "../00.data"
output_dir <- "../01.result/own_analysis"

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
## 14. GO CC Figure 5A-like plot
############################################################
## localization 파일과 별개로 GO Cellular Component 기반 plot도 생성
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

ont_use <- "CC"

go_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(gene_stat_go$ENTREZID),
  keytype = "ENTREZID",
  columns = c("GO", "ONTOLOGY")
) %>%
  filter(!is.na(GO), ONTOLOGY == ont_use) %>%
  distinct(ENTREZID, GO, .keep_all = TRUE)

go_term_names <- AnnotationDbi::select(
  GO.db,
  keys = unique(go_map$GO),
  keytype = "GOID",
  columns = c("TERM")
) %>%
  distinct(GOID, .keep_all = TRUE)

dat_go <- go_map %>%
  inner_join(gene_stat_go, by = "ENTREZID")

all_genes_go <- gene_stat_go %>%
  dplyr::select(ENTREZID, clip_enrichment, ribo_density_change) %>%
  distinct()

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
    clip_fdr = p.adjust(clip_p, method = "BH"),
    ribo_fdr = p.adjust(ribo_p, method = "BH"),
    combined_fdr = pmin(clip_fdr, ribo_fdr, na.rm = TRUE),
    neglog10_fdr = -log10(combined_fdr)
  ) %>%
  filter(
    n_genes >= 10,
    n_genes <= 3000,
    combined_fdr < 0.05
  ) %>%
  arrange(combined_fdr)

label_go <- go_stats %>%
  filter(
    combined_fdr < 1e-5 |
      str_detect(
        str_to_lower(TERM),
        "endoplasmic|membrane|golgi|extracellular|cytoplasm|nucleus|mitochond"
      )
  ) %>%
  slice_head(n = 30)

p_go <- ggplot(
  go_stats,
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
    alpha = 0.75
  ) +
  geom_text_repel(
    data = label_go,
    aes(label = paste0(TERM, " (", n_genes, ")")),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_size_continuous(
    range = c(2, 14),
    name = "Number of genes"
  ) +
  scale_color_gradient(
    low = "khaki1",
    high = "red4",
    name = "-log10(FDR)"
  ) +
  labs(
    title = "GO Cellular Component analysis for CLIP and ribosome profiling",
    x = "Mean LIN28A CLIP enrichment (log2)",
    y = "Mean ribosome density change upon Lin28a knockdown (log2)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_line(color = "grey85", linetype = "dashed"),
    plot.title = element_text(face = "bold")
  )

print(p_go)

ggsave(
  file.path(output_dir, "GO_CC_Figure5A_like_with_localization_context.png"),
  p_go,
  width = 12,
  height = 7,
  dpi = 900
)