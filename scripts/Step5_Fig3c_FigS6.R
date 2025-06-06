#!/usr/bin/env Rscript
library(parallel)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(gridExtra)
library(tidyverse)

mean_computing <- function(data) {
  data %>%
    group_by(V1) %>%
    mutate(mean_FC = mean(V4)) %>%
    ungroup() %>%
    na.omit()
}

#Set your actual path
setwd("")

#Hap1 comparison without window####
file_to_open <- list.files(pattern = "_AS_HOR_hap1.bed")
all_files_list <- mclapply(file_to_open, fread, mc.cores = 8)
names(all_files_list) <- file_to_open
all_files_list <- mclapply(all_files_list, mean_computing, mc.cores = 8)

all_files_list_df <- bind_rows(all_files_list)
all_files_list_df$identificator <- paste0(all_files_list_df$V1, "_", all_files_list_df$V5, "_", all_files_list_df$V6, "_", all_files_list_df$V7)
setDT(all_files_list_df)
plot_data <- all_files_list_df[, .SD[1], by = identificator]
plot_data$identificator2 <- paste0(plot_data$V5, "_", plot_data$V6)

level_order <- c(paste0("chr", 1:22, "_hap1"), "chrX_hap1")
plot_data$V1 <- factor(plot_data$V1, levels = level_order)

chromosome_colors <- c(
  chr1_hap1 = "#1b9e77",
  chr2_hap1 = "#d95f02",
  chr3_hap1 = "#7570b3",
  chr4_hap1 = "#e7298a",
  chr5_hap1 = "#66a61e",
  chr6_hap1 = "#e6ab02",
  chr7_hap1 = "#a6761d",
  chr8_hap1 = "#666666",
  chr9_hap1 = "#1f78b4",
  chr10_hap1= "#b2df8a",
  chr11_hap1 = "#fb9a99",
  chr12_hap1 = "#fdbf6f",
  chr13_hap1 = "#cab2d6",
  chr14_hap1 = "#ffff99",
  chr15_hap1 = "#a6cee3",
  chr16_hap1 = "#b15928",
  chr17_hap1 = "#8dd3c7",
  chr18_hap1 = "#ffffb3",
  chr19_hap1 = "#bebada",
  chr20_hap1 = "#fb8072",
  chr21_hap1 = "#80b1d3",
  chr22_hap1 = "#fdb462",
  chrX_hap1 = "#bc80bd"
)

#paired t-test for comparison across centromeres####
test_in <- plot_data[plot_data$V7 == "in_centromere", c("V1", "mean_FC", "V6", "V5", "identificator2")]
colnames(test_in) <- c("chromosome", "mean_FC", "condition", "replicate", "identificator2")
paired_data <- test_in[, -c("condition", "replicate")] %>% tidyr::pivot_wider(names_from = identificator2, values_from = mean_FC)
sum(!complete.cases(paired_data))
paired_data <- na.omit(paired_data)

t_test_results_in_hap1 <- paired_data %>%
  group_by(chromosome) %>%
  summarise(
    t_test = list(
      t.test(
        x = c(REP1_siEP, REP2_siEP, REP3_siEP),
        y = c(REP1_siNG, REP2_siNG, REP3_siNG),
        paired = F
      )
    ),
    # Compute effect size: Cohen's d for paired samples
    cohens_d = {
      diffs <- c(REP1_siEP - REP1_siNG, REP2_siEP - REP2_siNG, REP3_siEP - REP3_siNG)
      mean(diffs) / sd(diffs)
    }
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    mean_diff = sapply(t_test, function(x) x$estimate[1] - x$estimate[2]),  #for paired: sapply(t_test, function(x) x$estimate)
    statistic = sapply(t_test, function(x) x$statistic),
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_pval = -log10(p_value)
  ) %>%
  select(chromosome, p_value, p_adj, neg_log10_pval, mean_diff, statistic, cohens_d)

# Save your summarized data of in_centromere signals
summary_data <- plot_data %>%
  filter(V7 == "in_centromere") %>%
  group_by(V1, V6, V5) %>%
  summarise(value = mean_FC, .groups = "drop") %>%
  group_by(V1, V6) %>%
  summarise(
    mean_FC = mean(value),
    se_FC = sd(value) / sqrt(n()),
    .groups = "drop"
  )

level_order_condition <- c("siNG", "siEP")
summary_data$V6 <- factor(summary_data$V6, levels = level_order_condition)

#Position of the asterisk
asterisk_data <- summary_data %>%
  group_by(V1) %>%
  summarise(
    y_pos = max(mean_FC + se_FC) + 0.02
  )

#Table for asterisk
asterisk_data <- data.frame(
  V1 = c(paste0("chr", 1:22, "_hap1"), "chrX_hap1"),
  y_pos = asterisk_data$y_pos,
  label = c("*", "**", "**", "**", "*", "*", "*", "**", "*", "**", "*", "**", "**", "**", "**", "**", "**", "*", "*", "*", "**", "*", "**")
)

pdf("by_chr_in_centromeres_hap1.pdf", width = 20, height = 8)
# Plot
ggplot(summary_data, aes(x = V1, y = mean_FC, color = V6)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 1.5
  ) +
  geom_errorbar(
    aes(ymin = mean_FC - se_FC, ymax = mean_FC + se_FC),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = asterisk_data,
    aes(x = V1, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  scale_color_manual(
    values = c("siEP"= "#ff3333", "siNG" = "#3333ff"),
    name = "Condition"
  ) +
  labs(
    x = "",
    y = "mean FC ± SE",
    title = "CENP-A enrichment within centromeres Hap1 (mean of FC in the three replicates per condition ± SE)"
  ) +
  theme_linedraw(base_size = 9) +
  ylim(0,1.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
dev.off()

#paired t-test for comparison outside of centromeres####
test_in <- plot_data[plot_data$V7 == "out_centromere", c("V1", "mean_FC", "V6", "V5", "identificator2")]
colnames(test_in) <- c("chromosome", "mean_FC", "condition", "replicate", "identificator2")
paired_data <- test_in[, -c("condition", "replicate")] %>% tidyr::pivot_wider(names_from = identificator2, values_from = mean_FC)
sum(!complete.cases(paired_data))
paired_data <- na.omit(paired_data)

t_test_results_out_hap1 <- paired_data %>%
  group_by(chromosome) %>%
  summarise(
    t_test = list(
      t.test(
        x = c(REP1_siEP, REP2_siEP, REP3_siEP),
        y = c(REP1_siNG, REP2_siNG, REP3_siNG),
        paired = F
      )
    ),
    # Compute effect size: Cohen's d for paired samples
    cohens_d = {
      diffs <- c(REP1_siEP - REP1_siNG, REP2_siEP - REP2_siNG, REP3_siEP - REP3_siNG)
      mean(diffs) / sd(diffs)
    }
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    mean_diff = sapply(t_test, function(x) x$estimate[1] - x$estimate[2]),  #for paired: sapply(t_test, function(x) x$estimate)
    statistic = sapply(t_test, function(x) x$statistic),
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_pval = -log10(p_value)
  ) %>%
  select(chromosome, p_value, p_adj, neg_log10_pval, mean_diff, statistic, cohens_d)

# Save your summarized data of out_centromere signals
summary_data <- plot_data %>%
  filter(V7 == "out_centromere") %>%
  group_by(V1, V6, V5) %>%
  summarise(value = mean_FC, .groups = "drop") %>%
  group_by(V1, V6) %>%
  summarise(
    mean_FC = mean(value),
    se_FC = sd(value) / sqrt(n()),
    .groups = "drop"
  )

level_order_condition <- c("siNG", "siEP")
summary_data$V6 <- factor(summary_data$V6, levels = level_order_condition)

asterisk_data <- summary_data %>%
  group_by(V1) %>%
  summarise(
    y_pos = max(mean_FC + se_FC) + 0.02  # Adjust offset as needed
  )

pdf("by_chr_out_centromeres_hap1.pdf", width = 20, height = 8)
# Plot
ggplot(summary_data, aes(x = V1, y = mean_FC, color = V6)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 1.5
  ) +
  geom_errorbar(
    aes(ymin = mean_FC - se_FC, ymax = mean_FC + se_FC),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = asterisk_data,
    aes(x = V1, y = y_pos),
    label = "*",
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  scale_color_manual(
    values = c("siEP"= "#ff3333", "siNG" = "#3333ff"),
    name = "Condition"
  ) +
  labs(
    x = "",
    y = "mean FC ± SE",
    title = "CENP-A enrichment outside centromeres Hap1 (mean of FC in the three replicates per condition ± SE)"
  ) +
  theme_linedraw(base_size = 9) +
  ylim(0,1.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
dev.off()

#mean of the mean####
data_mean <- plot_data
data_mean$identificator3 <- paste0(data_mean$V5, "_", data_mean$V6, "_", data_mean$V7)
data_mean <- data_mean %>% group_by(identificator3) %>% mutate(total_mean_FC = mean(mean_FC)) %>% ungroup()
setDT(data_mean)
data_mean <- data_mean[, .SD[1], by = identificator3]

data_mean_f <- data_mean[, c("V7", "total_mean_FC", "identificator2")]
data_mean_f <- as.data.frame(data_mean_f)

data_mean_f <- data_mean_f %>% tidyr::pivot_wider(names_from = identificator2, values_from = total_mean_FC)

t_test_results_hap1_total <- data_mean_f %>%
  group_by(V7) %>%
  summarise(
    t_test = list(
      t.test(
        x = c(REP1_siEP, REP2_siEP, REP3_siEP),
        y = c(REP1_siNG, REP2_siNG, REP3_siNG),
        paired = F
      )
    ),
    # Compute effect size: Cohen's d for paired samples
    cohens_d = {
      diffs <- c(REP1_siEP - REP1_siNG, REP2_siEP - REP2_siNG, REP3_siEP - REP3_siNG)
      mean(diffs) / sd(diffs)
    }
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    mean_diff = sapply(t_test, function(x) x$estimate[1] - x$estimate[2]),  #for paired: sapply(t_test, function(x) x$estimate)
    statistic = sapply(t_test, function(x) x$statistic),
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_pval = -log10(p_value)
  ) %>%
  select(V7, p_value, p_adj, neg_log10_pval, mean_diff, statistic, cohens_d)

summary_data <- data_mean %>%
  group_by(V7, V6, V5) %>%
  summarise(value = total_mean_FC, .groups = "drop") %>%
  group_by(V7, V6) %>%
  summarise(
    mean_FC = mean(value),
    se_FC = sd(value) / sqrt(n()),
    .groups = "drop"
  )

level_order_condition <- c("siNG", "siEP")
summary_data$V6 <- factor(summary_data$V6, levels = level_order_condition)

asterisk_data <- summary_data %>%
  group_by(V7) %>%
  summarise(
    y_pos = max(mean_FC + se_FC) + 0.02  # Adjust offset as needed
  )

# Creo i dati per gli asterischi
asterisk_data <- data.frame(
  V7 = c("out_centromere", "in_centromere"),
  y_pos = c(0.404, 0.778),   # posizioni verticali (devi adattare in base al tuo plot)
  label = c("*", "**")    # 1 asterisco per out, 2 per in
)

pdf("total_mean_centr_vs_nocentr_hap1.pdf", width = 20, height = 8)
# Plot
ggplot(summary_data, aes(x = V7, y = mean_FC, color = V6)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 1.5
  ) +
  geom_errorbar(
    aes(ymin = mean_FC - se_FC, ymax = mean_FC + se_FC),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = asterisk_data,
    aes(x = V7, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  scale_color_manual(
    values = c("siEP"= "#ff3333", "siNG" = "#3333ff"),
    name = "Condition"
  ) +
  labs(
    x = "",
    y = "mean FC ± SE",
    title = "Genome-wide CENP-A enrichment in Hap1"
  ) +
  theme_linedraw(base_size = 9) +
  ylim(0,1.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 8)
  )
dev.off()

#Hap2 comparison without window####
mean_computing <- function(data) {
  data %>%
    group_by(V1) %>%
    mutate(mean_FC = mean(V4)) %>%
    ungroup() %>%
    na.omit()
}

file_to_open <- list.files(pattern = "_AS_HOR_hap2.bed")
all_files_list <- mclapply(file_to_open, fread, mc.cores = 8)
names(all_files_list) <- file_to_open
all_files_list <- mclapply(all_files_list, mean_computing, mc.cores = 8)

all_files_list_df <- bind_rows(all_files_list)
all_files_list_df$identificator <- paste0(all_files_list_df$V1, "_", all_files_list_df$V5, "_", all_files_list_df$V6, "_", all_files_list_df$V7)
setDT(all_files_list_df)
plot_data <- all_files_list_df[, .SD[1], by = identificator]
plot_data$identificator2 <- paste0(plot_data$V5, "_", plot_data$V6)

level_order <- c(paste0("chr", 1:22, "_hap2"), "chrX_hap2")
plot_data$V1 <- factor(plot_data$V1, levels = level_order)

chromosome_colors <- c(
  chr1_hap2 = "#1b9e77",
  chr2_hap2 = "#d95f02",
  chr3_hap2 = "#7570b3",
  chr4_hap2 = "#e7298a",
  chr5_hap2 = "#66a61e",
  chr6_hap2 = "#e6ab02",
  chr7_hap2 = "#a6761d",
  chr8_hap2 = "#666666",
  chr9_hap2 = "#1f78b4",
  chr10_hap2= "#b2df8a",
  chr11_hap2 = "#fb9a99",
  chr12_hap2 = "#fdbf6f",
  chr13_hap2 = "#cab2d6",
  chr14_hap2 = "#ffff99",
  chr15_hap2 = "#a6cee3",
  chr16_hap2 = "#b15928",
  chr17_hap2 = "#8dd3c7",
  chr18_hap2 = "#ffffb3",
  chr19_hap2 = "#bebada",
  chr20_hap2 = "#fb8072",
  chr21_hap2 = "#80b1d3",
  chr22_hap2 = "#fdb462",
  chrX_hap2 = "#bc80bd"
)

#paired t-test for comparison across centromeres####
test_in <- plot_data[plot_data$V7 == "in_centromere", c("V1", "mean_FC", "V6", "V5", "identificator2")]
colnames(test_in) <- c("chromosome", "mean_FC", "condition", "replicate", "identificator2")
paired_data <- test_in[, -c("condition", "replicate")] %>% tidyr::pivot_wider(names_from = identificator2, values_from = mean_FC)
sum(!complete.cases(paired_data))
paired_data <- na.omit(paired_data)

t_test_results_in_hap2 <- paired_data %>%
  group_by(chromosome) %>%
  summarise(
    t_test = list(
      t.test(
        x = c(REP1_siEP, REP2_siEP, REP3_siEP),
        y = c(REP1_siNG, REP2_siNG, REP3_siNG),
        paired = F
      )
    ),
    # Compute effect size: Cohen's d for paired samples
    cohens_d = {
      diffs <- c(REP1_siEP - REP1_siNG, REP2_siEP - REP2_siNG, REP3_siEP - REP3_siNG)
      mean(diffs) / sd(diffs)
    }
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    mean_diff = sapply(t_test, function(x) x$estimate[1] - x$estimate[2]),  #for paired: sapply(t_test, function(x) x$estimate)
    statistic = sapply(t_test, function(x) x$statistic),
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_pval = -log10(p_value)
  ) %>%
  select(chromosome, p_value, p_adj, neg_log10_pval, mean_diff, statistic, cohens_d)

# Save your summarized data
summary_data <- plot_data %>%
  filter(V7 == "in_centromere") %>%
  group_by(V1, V6, V5) %>%
  summarise(value = mean_FC, .groups = "drop") %>%
  group_by(V1, V6) %>%
  summarise(
    mean_FC = mean(value),
    se_FC = sd(value) / sqrt(n()),
    .groups = "drop"
  )

level_order_condition <- c("siNG", "siEP")
summary_data$V6 <- factor(summary_data$V6, levels = level_order_condition)

asterisk_data <- summary_data %>%
  group_by(V1) %>%
  summarise(
    y_pos = max(mean_FC + se_FC) + 0.02  # Adjust offset as needed
  )

# Creo i dati per gli asterischi
asterisk_data <- data.frame(
  V1 = c(paste0("chr", 1:22, "_hap2"), "chrX_hap2"),
  y_pos = asterisk_data$y_pos,   # posizioni verticali (devi adattare in base al tuo plot)
  label = c("*", "**", "**", "**", "**", "*", "*", "**", "*", "**", "*", "**", "**", "**", "**", "**", "*", "*", "*", "*", "*", "*", "**")    # 1 asterisco per out, 2 per in
)

pdf("by_chr_in_centromeres_hap2.pdf", width = 20, height = 8)
# Plot
ggplot(summary_data, aes(x = V1, y = mean_FC, color = V6)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 1.5
  ) +
  geom_errorbar(
    aes(ymin = mean_FC - se_FC, ymax = mean_FC + se_FC),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = asterisk_data,
    aes(x = V1, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  scale_color_manual(
    values = c("siEP"= "#ff3333", "siNG" = "#3333ff"),
    name = "Condition"
  ) +
  labs(
    x = "",
    y = "mean FC ± SE",
    title = "CENP-A enrichment within centromeres Hap2 (mean of FC in the three replicates per condition ± SE)"
  ) +
  theme_linedraw(base_size = 9) +
  ylim(0,1.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
dev.off()

#paired t-test for comparison outside of centromeres####
test_in <- plot_data[plot_data$V7 == "out_centromere", c("V1", "mean_FC", "V6", "V5", "identificator2")]
colnames(test_in) <- c("chromosome", "mean_FC", "condition", "replicate", "identificator2")
paired_data <- test_in[, -c("condition", "replicate")] %>% tidyr::pivot_wider(names_from = identificator2, values_from = mean_FC)
sum(!complete.cases(paired_data))
paired_data <- na.omit(paired_data)

t_test_results_out_hap2 <- paired_data %>%
  group_by(chromosome) %>%
  summarise(
    t_test = list(
      t.test(
        x = c(REP1_siEP, REP2_siEP, REP3_siEP),
        y = c(REP1_siNG, REP2_siNG, REP3_siNG),
        paired = F
      )
    ),
    # Compute effect size: Cohen's d for paired samples
    cohens_d = {
      diffs <- c(REP1_siEP - REP1_siNG, REP2_siEP - REP2_siNG, REP3_siEP - REP3_siNG)
      mean(diffs) / sd(diffs)
    }
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    mean_diff = sapply(t_test, function(x) x$estimate[1] - x$estimate[2]),  #for paired: sapply(t_test, function(x) x$estimate)
    statistic = sapply(t_test, function(x) x$statistic),
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_pval = -log10(p_value)
  ) %>%
  select(chromosome, p_value, p_adj, neg_log10_pval, mean_diff, statistic, cohens_d)

# Save your summarized data
summary_data <- plot_data %>%
  filter(V7 == "out_centromere") %>%
  group_by(V1, V6, V5) %>%
  summarise(value = mean_FC, .groups = "drop") %>%
  group_by(V1, V6) %>%
  summarise(
    mean_FC = mean(value),
    se_FC = sd(value) / sqrt(n()),
    .groups = "drop"
  )

level_order_condition <- c("siNG", "siEP")
summary_data$V6 <- factor(summary_data$V6, levels = level_order_condition)

asterisk_data <- summary_data %>%
  group_by(V1) %>%
  summarise(
    y_pos = max(mean_FC + se_FC) + 0.02  # Adjust offset as needed
  )

pdf("by_chr_out_centromeres_hap2.pdf", width = 20, height = 8)
# Plot
ggplot(summary_data, aes(x = V1, y = mean_FC, color = V6)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 1.5
  ) +
  geom_errorbar(
    aes(ymin = mean_FC - se_FC, ymax = mean_FC + se_FC),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = asterisk_data,
    aes(x = V1, y = y_pos),
    label = "*",
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  scale_color_manual(
    values = c("siEP"= "#ff3333", "siNG" = "#3333ff"),
    name = "Condition"
  ) +
  labs(
    x = "",
    y = "mean FC ± SE",
    title = "CENP-A enrichment outside centromeres Hap2 (mean of FC in the three replicates per condition ± SE)"
  ) +
  theme_linedraw(base_size = 9) +
  ylim(0,1.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
dev.off()

#mean of the mean####
data_mean <- plot_data
data_mean$identificator3 <- paste0(data_mean$V5, "_", data_mean$V6, "_", data_mean$V7)
data_mean <- data_mean %>% group_by(identificator3) %>% mutate(total_mean_FC = mean(mean_FC)) %>% ungroup()
setDT(data_mean)
data_mean <- data_mean[, .SD[1], by = identificator3]

data_mean_f <- data_mean[, c("V7", "total_mean_FC", "identificator2")]
data_mean_f <- as.data.frame(data_mean_f)

data_mean_f <- data_mean_f %>% tidyr::pivot_wider(names_from = identificator2, values_from = total_mean_FC)

t_test_results_hap2_total <- data_mean_f %>%
  group_by(V7) %>%
  summarise(
    t_test = list(
      t.test(
        x = c(REP1_siEP, REP2_siEP, REP3_siEP),
        y = c(REP1_siNG, REP2_siNG, REP3_siNG),
        paired = F
      )
    ),
    # Compute effect size: Cohen's d for paired samples
    cohens_d = {
      diffs <- c(REP1_siEP - REP1_siNG, REP2_siEP - REP2_siNG, REP3_siEP - REP3_siNG)
      mean(diffs) / sd(diffs)
    }
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    mean_diff = sapply(t_test, function(x) x$estimate[1] - x$estimate[2]),  #for paired: sapply(t_test, function(x) x$estimate)
    statistic = sapply(t_test, function(x) x$statistic),
    p_adj = p.adjust(p_value, method = "BH"),
    neg_log10_pval = -log10(p_value)
  ) %>%
  select(V7, p_value, p_adj, neg_log10_pval, mean_diff, statistic, cohens_d)

summary_data <- data_mean %>%
  group_by(V7, V6, V5) %>%
  summarise(value = total_mean_FC, .groups = "drop") %>%
  group_by(V7, V6) %>%
  summarise(
    mean_FC = mean(value),
    se_FC = sd(value) / sqrt(n()),
    .groups = "drop"
  )

level_order_condition <- c("siNG", "siEP")
summary_data$V6 <- factor(summary_data$V6, levels = level_order_condition)

asterisk_data <- summary_data %>%
  group_by(V7) %>%
  summarise(
    y_pos = max(mean_FC + se_FC) + 0.02  # Adjust offset as needed
  )

# Creo i dati per gli asterischi
asterisk_data <- data.frame(
  V7 = c("out_centromere", "in_centromere"),
  y_pos = c(0.404, 0.778),   # posizioni verticali (devi adattare in base al tuo plot)
  label = c("*", "**")    # 1 asterisco per out, 2 per in
)

pdf("total_mean_centr_vs_nocentr_hap2.pdf", width = 20, height = 8)
# Plot
ggplot(summary_data, aes(x = V7, y = mean_FC, color = V6)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 1.5
  ) +
  geom_errorbar(
    aes(ymin = mean_FC - se_FC, ymax = mean_FC + se_FC),
    width = 0.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = asterisk_data,
    aes(x = V7, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  scale_color_manual(
    values = c("siEP"= "#ff3333", "siNG" = "#3333ff"),
    name = "Condition"
  ) +
  labs(
    x = "",
    y = "mean FC ± SE",
    title = "Genome-wide CENP-A enrichment in Hap2"
  ) +
  theme_linedraw(base_size = 9) +
  ylim(0,1.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 8)
  )
dev.off()