library(tidyverse)
library(dplyr)
library(ggplot2)

resultpath <- "~/GuidedMission/01.result"
setwd(resultpath)

read_start_file <- function(file, sample_name) {
  read_tsv(
    file,
    col_names = c(
      "chr",
      "fivep_start",
      "fivep_end",
      "count",
      "exon_chr",
      "exon_start",
      "exon_end",
      "transcript_id",
      "codon_start",
      "strand"
    ),
    show_col_types = FALSE
  ) %>%
    mutate(
      sample = sample_name,
      rel_pos = fivep_start - codon_start
    ) %>%
    filter(rel_pos >= -50, rel_pos <= 50) %>%
    group_by(sample, rel_pos) %>%
    summarise(read_count = sum(count), .groups = "drop") %>%
    mutate(read_count_x1000 = read_count / 1000)
}

start_df <- bind_rows(
  read_start_file("fivepcounts-filtered-RPF-siLuc.txt", "siLuc"),
  read_start_file("fivepcounts-filtered-RPF-siLin28a.txt", "siLin28a")
)

p_start <- ggplot(start_df, aes(x = rel_pos, y = read_count_x1000)) +
  geom_col(width = 0.8, fill = "grey20") +
  geom_vline(xintercept = 0, color = "firebrick", linewidth = 0.7) +
  facet_grid(sample ~ ., scales = "fixed") +
  scale_x_continuous(
    limits = c(-50, 50),
    breaks = seq(-50, 50, 10)
  ) +
  scale_y_continuous(
    limits = c(0, 120),
    breaks = seq(0, 120, 40),
    minor_breaks = seq(0, 120, 20),
    expand = c(0, 0)
  ) +
  labs(
    title = "Ribosome footprint density near start codons",
    x = "Relative position to start codon of 5′-end of reads",
    y = "Raw read count\n(x1000)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.text.y = element_text(size = 18, face = "bold", angle = 0),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor.y = element_line(color = "grey85"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 20, face = "bold")
  )

ggsave(
  filename = "RPF_density_near_start_codon.png",
  plot = p_start,
  width = 8,
  height = 5,
  dpi = 300
)