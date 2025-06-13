library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)

setwd("/data/VM/Documents/mylab/CENPA_Basrai_Munira/F7796_vx0054_2/BAM_SPIKEIN/BigWig_SpikeIn/intersect_in_out_centromere/prova_github_codes/Fig3d_FigS7b")

plot_log2_conditions <- function(title) {
  ggplot(mean_data, aes(x = V1, y = mean_value, color = regions)) +
    geom_point(position = position_dodge(width = 0.5), size = 1.5, alpha = 0.8) +
    ggtitle(title) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    ylab("Mean LogFC(siEP/siNG)") +
    xlab("") +
    ylim(-0.5, 1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),  # Rotate x-axis labels by 45 degrees
      axis.text.y = element_text(colour = "black"),
      axis.title.y = element_text(size = 8),# Y-axis labels in black
      plot.margin = margin(10, 10, 10, 10)  # Adjust plot margins to make it more compact
    ) +
    scale_color_manual(values = c(
      "Live" = "red4",        # Customize the color for 'Live'
      "Divergent" = "yellow4",    # Customize the color for 'Divergent'
      "CDRs" = "green4",       # Customize the color for 'CDRs'
      "Inactive" = "orange4"   # Customize the color for 'Inactive'
    )) +
    coord_cartesian(clip = "off") + theme(legend.position = "none")
}

#HAP1_REP1####
FC_by_regions <- fread("log2_rep1_hap1.bedgraph_AS-HOR.bed", header = F)
FC_by_regions <- FC_by_regions %>%
  mutate(regions = case_when(
    grepl("L", V4) ~ "Live",
    grepl("^cdr_\\d", V4) ~ "CDRs",
    grepl("d", V4) ~ "Divergent",
    TRUE ~ "Inactive"
  ))

mean_data <- FC_by_regions %>%
  group_by(V1, regions) %>%
  summarise(mean_value = mean(V8), .groups = 'drop')

mean_data <- mean_data %>%
  mutate(V1 = factor(V1, levels = paste0("chr", c(1:22, "X"), "_hap1")))


p1 <- plot_log2_conditions("Hap1 Rep1")

#HAP1_REP2####
FC_by_regions <- fread("log2_rep2_hap1.bedgraph_AS-HOR.bed", header = F)
FC_by_regions <- FC_by_regions %>%
  mutate(regions = case_when(
    grepl("L", V4) ~ "Live",
    grepl("^cdr_\\d", V4) ~ "CDRs",
    grepl("d", V4) ~ "Divergent",
    TRUE ~ "Inactive"
  ))

mean_data <- FC_by_regions %>%
  group_by(V1, regions) %>%
  summarise(mean_value = mean(V8), .groups = 'drop')

mean_data <- mean_data %>%
  mutate(V1 = factor(V1, levels = paste0("chr", c(1:22, "X"), "_hap1")))


p2 <- plot_log2_conditions("Hap1 Rep2")

#HAP1_REP3####
FC_by_regions <- fread("log2_rep3_hap1.bedgraph_AS-HOR.bed", header = F)
FC_by_regions <- FC_by_regions %>%
  mutate(regions = case_when(
    grepl("L", V4) ~ "Live",
    grepl("^cdr_\\d", V4) ~ "CDRs",
    grepl("d", V4) ~ "Divergent",
    TRUE ~ "Inactive"
  ))

mean_data <- FC_by_regions %>%
  group_by(V1, regions) %>%
  summarise(mean_value = mean(V8), .groups = 'drop')

mean_data <- mean_data %>%
  mutate(V1 = factor(V1, levels = paste0("chr", c(1:22, "X"), "_hap1")))


p3 <- plot_log2_conditions("Hap1 Rep3")

#HAP2_REP1####
FC_by_regions <- fread("log2_rep1_hap2.bedgraph_AS-HOR.bed", header = F)
FC_by_regions <- FC_by_regions %>%
  mutate(regions = case_when(
    grepl("L", V4) ~ "Live",
    grepl("^cdr_\\d", V4) ~ "CDRs",
    grepl("d", V4) ~ "Divergent",
    TRUE ~ "Inactive"
  ))

mean_data <- FC_by_regions %>%
  group_by(V1, regions) %>%
  summarise(mean_value = mean(V8), .groups = 'drop')

mean_data <- mean_data %>%
  mutate(V1 = factor(V1, levels = paste0("chr", c(1:22, "X"), "_hap2")))


p4 <- plot_log2_conditions("Hap2 Rep1")

#HAP2_REP2####
FC_by_regions <- fread("log2_rep2_hap2.bedgraph_AS-HOR.bed", header = F)
FC_by_regions <- FC_by_regions %>%
  mutate(regions = case_when(
    grepl("L", V4) ~ "Live",
    grepl("^cdr_\\d", V4) ~ "CDRs",
    grepl("d", V4) ~ "Divergent",
    TRUE ~ "Inactive"
  ))

mean_data <- FC_by_regions %>%
  group_by(V1, regions) %>%
  summarise(mean_value = mean(V8), .groups = 'drop')

mean_data <- mean_data %>%
  mutate(V1 = factor(V1, levels = paste0("chr", c(1:22, "X"), "_hap2")))


p5 <- plot_log2_conditions("Hap2 Rep2")

#HAP2_REP3####
FC_by_regions <- fread("log2_rep3_hap2.bedgraph_AS-HOR.bed", header = F)
FC_by_regions <- FC_by_regions %>%
  mutate(regions = case_when(
    grepl("L", V4) ~ "Live",
    grepl("^cdr_\\d", V4) ~ "CDRs",
    grepl("d", V4) ~ "Divergent",
    TRUE ~ "Inactive"
  ))

mean_data <- FC_by_regions %>%
  group_by(V1, regions) %>%
  summarise(mean_value = mean(V8), .groups = 'drop')

mean_data <- mean_data %>%
  mutate(V1 = factor(V1, levels = paste0("chr", c(1:22, "X"), "_hap2")))


p6 <- plot_log2_conditions("Hap2 Rep3")

pdf("log2_by_replicate.pdf", width = 20, height = 8)
grid.arrange(p1,p4,p2,p5,p3,p6, ncol = 2)
dev.off()