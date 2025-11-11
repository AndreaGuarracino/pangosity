#!/usr/bin/env Rscript

library(tidyverse)
library(patchwork)

# Disable scientific notation
options(scipen = 999)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

# All arguments are now required
if (length(args) < 6) {
  cat("Usage: Rscript plot_pangosity.R <samples.txt> <zygosity.tsv> <node_lengths.tsv> <output.pdf> <panplexity_mask.txt> <node_coverage_mask.txt>\n")
  cat("  samples.txt: Input file with sample_name<tab>coverage_file\n")
  cat("  zygosity.tsv: Output zygosity matrix from pangosity\n")
  cat("  node_lengths.tsv: 2-column file with node_id<tab>length\n")
  cat("  output.pdf: Output plot file\n")
  cat("  panplexity_mask.txt: File with mask values (0=low-complexity, 1=not masked), one per line per node\n")
  cat("  node_coverage_mask.txt: File with node coverage mask values, one per line per node\n")
  quit(status = 1)
}

samples_file <- args[1]
zygosity_file <- args[2]
node_lengths_file <- args[3]
output_pdf <- args[4]
panplexity_mask_file <- args[5]
node_coverage_mask_file <- args[6]

cat("Reading zygosity matrix:", zygosity_file, "\n")
# Read zygosity matrix
zyg <- read_tsv(zygosity_file, comment = "", show_col_types = FALSE, col_types = cols(.default = "c"))
colnames(zyg)[1] <- "node"
# Convert node column to numeric
zyg$node <- as.numeric(zyg$node)

cat("Reading sample list and coverage files\n")
# Read sample list and coverage
samples <- read_tsv(samples_file, col_names = c("sample", "coverage_file"), show_col_types = FALSE)

coverage_list <- lapply(1:nrow(samples), function(i) {
  cat("  Processing:", samples$coverage_file[i], "\n")
  lines <- readLines(samples$coverage_file[i])
  
  # Detect file format
  if (length(lines) > 0 && grepl("^##sample:", lines[1])) {
    # Pack format: ##sample: name\n#coverage\nvalue1\nvalue2\n...
    cat("    Detected pack format\n")
    
    # Filter out header lines (starting with #) and empty lines
    data_lines <- lines[!grepl("^#", lines) & lines != ""]
    # Changed from 0:(length(data_lines)-1) to 1:length(data_lines)
    tibble(sample = samples$sample[i], node = 1:length(data_lines), coverage = as.numeric(data_lines))
    
  } else if (length(lines) > 0 && grepl("^#sample\\t", lines[1])) {
    # Row format: #sample\tnode.1\tnode.2\t...
    cat("    Detected row format\n")
    cov_data <- read_tsv(samples$coverage_file[i], comment = "", show_col_types = FALSE)
    sample_name <- cov_data[[1]][1]
    cov_values <- as.numeric(cov_data[1, -1])
    # Changed from 0:(length(cov_values)-1) to 1:length(cov_values)
    tibble(sample = samples$sample[i], node = 1:length(cov_values), coverage = cov_values)
    
  } else {
    # Column format: one value per line (with optional # comments)
    cat("    Detected column format\n")
    data_lines <- lines[!grepl("^#", lines) & lines != ""]
    # Changed from 0:(length(data_lines)-1) to 1:length(data_lines)
    tibble(sample = samples$sample[i], node = 1:length(data_lines), coverage = as.numeric(data_lines))
  }
})
coverage <- bind_rows(coverage_list)

# Read node lengths (now required)
cat("Reading node lengths:", node_lengths_file, "\n")
node_lens <- read_tsv(node_lengths_file, col_names = c("node", "length"), show_col_types = FALSE)
node_lens <- node_lens %>% 
  arrange(node) %>%
  mutate(
    cum_start = cumsum(lag(length, default = 0)),
    cum_end = cum_start + length
  )

# Read panplexity mask (now required)
cat("Reading panplexity mask:", panplexity_mask_file, "\n")
mask_lines <- readLines(panplexity_mask_file)
mask_lines <- mask_lines[!grepl("^#", mask_lines) & mask_lines != ""]
mask_values <- as.integer(mask_lines)

panplexity_mask <- tibble(
  node = 1:length(mask_values),  # Changed from 0:(length(mask_values)-1)
  masked = mask_values
) %>%
  mutate(masked_label = ifelse(masked == 0, "Low-complexity", "Not masked")) %>%
  left_join(node_lens %>% select(node, cum_start, cum_end, length), by = "node")

cat("  Low-complexity nodes:", sum(mask_values == 0), "/", length(mask_values), "\n")
cat("  Not masked nodes:", sum(mask_values == 1), "/", length(mask_values), "\n")

# Read node coverage mask (new required parameter)
cat("Reading node coverage mask:", node_coverage_mask_file, "\n")
coverage_mask_lines <- readLines(node_coverage_mask_file)
coverage_mask_lines <- coverage_mask_lines[!grepl("^#", coverage_mask_lines) & coverage_mask_lines != ""]
coverage_mask_values <- as.integer(coverage_mask_lines)

node_coverage_mask <- tibble(
  node = 1:length(coverage_mask_values),  # Changed from 0:(length(coverage_mask_values)-1)
  coverage_masked = coverage_mask_values
) %>%
  mutate(coverage_masked_label = ifelse(coverage_masked == 0, "Masked", "Not masked")) %>%
  left_join(node_lens %>% select(node, cum_start, cum_end, length), by = "node")

cat("  Masked nodes:", sum(coverage_mask_values == 0), "/", length(coverage_mask_values), "\n")
cat("  Not masked nodes:", sum(coverage_mask_values == 1), "/", length(coverage_mask_values), "\n")

# Join node positions
coverage <- coverage %>% 
  left_join(node_lens %>% select(node, cum_start, cum_end, length), by = "node")

zyg_long <- zyg %>% 
  pivot_longer(-node, names_to = "sample", values_to = "genotype") %>%
  left_join(node_lens %>% select(node, cum_start, cum_end, length), by = "node")

# Convert genotypes to numeric for visualization
zyg_long <- zyg_long %>%
  mutate(gt_num = case_when(
    genotype %in% c("0", "0/0") ~ 0,
    genotype %in% c("0/1") ~ 1,
    genotype %in% c("1", "1/1") ~ 2,
    TRUE ~ NA_real_
  ))

# Compute statistics
cat("Computing statistics\n")
zyg_stats <- zyg_long %>% 
  filter(!is.na(genotype)) %>%
  mutate(genotype_label = case_when(
    genotype == "0" ~ "0 (absent)",
    genotype == "1" ~ "1 (present)",
    genotype == "0/0" ~ "0/0 (hom ref)",
    genotype == "0/1" ~ "0/1 (het)",
    genotype == "1/1" ~ "1/1 (hom alt)",
    genotype %in% c(".", "./.") ~ "Missing",
    TRUE ~ as.character(genotype)
  )) %>%
  group_by(sample, genotype_label) %>%
  summarise(count = n(), .groups = "drop")

zyg_stats_pct <- zyg_stats %>%
  group_by(sample) %>%
  mutate(
    total = sum(count),
    percentage = count / total * 100,
    label = ifelse(percentage > 0, sprintf("%.1f%%", percentage), "")
  )

# Compute statistics by sequence length
cat("Computing sequence length statistics\n")
zyg_stats_len <- zyg_long %>% 
  filter(!is.na(genotype)) %>%
  mutate(genotype_label = case_when(
    genotype == "0" ~ "0 (absent)",
    genotype == "1" ~ "1 (present)",
    genotype == "0/0" ~ "0/0 (hom ref)",
    genotype == "0/1" ~ "0/1 (het)",
    genotype == "1/1" ~ "1/1 (hom alt)",
    genotype %in% c(".", "./.") ~ "Missing",
    TRUE ~ as.character(genotype)
  )) %>%
  group_by(sample, genotype_label) %>%
  summarise(total_length = sum(length, na.rm = TRUE), .groups = "drop")

zyg_stats_len_pct <- zyg_stats_len %>%
  group_by(sample) %>%
  mutate(
    total = sum(total_length),
    percentage = total_length / total * 100,
    label = ifelse(percentage > 0, sprintf("%.1f%%", percentage), "")
  )

cov_stats <- coverage %>%
  group_by(sample) %>%
  summarise(
    mean_cov = mean(coverage, na.rm = TRUE), 
    median_cov = median(coverage, na.rm = TRUE),
    sd_cov = sd(coverage, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate the global x-axis limits for all node-based plots
x_min <- 0
x_max <- max(node_lens$cum_end, na.rm = TRUE)

cat("Creating plots\n")

# Plot 1: Genotype Distribution by Sample (stacked bar chart)
p1 <- ggplot(zyg_stats_pct, aes(x = sample, y = count, fill = genotype_label)) +
  geom_col(position = "stack") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(
    values = c(
      "0 (absent)" = "#2166ac", "1 (present)" = "#b2182b",
      "0/0 (hom ref)" = "#2166ac", "0/1 (het)" = "#fddbc7", "1/1 (hom alt)" = "#b2182b",
      "Missing" = "grey80"),
    name = "Genotype",
    guide = guide_legend(override.aes = list(color = "black"))
  ) +
  labs(title = "Genotype distribution by sample (node count)", x = NULL, y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

# Plot 1b: Genotype Distribution by Sequence Length
p1b <- ggplot(zyg_stats_len_pct, aes(x = sample, y = total_length, fill = genotype_label)) +
  geom_col(position = "stack") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(
    values = c(
      "0 (absent)" = "#2166ac", "1 (present)" = "#b2182b",
      "0/0 (hom ref)" = "#2166ac", "0/1 (het)" = "#fddbc7", "1/1 (hom alt)" = "#b2182b",
      "Missing" = "grey80"),
    name = "Genotype",
    guide = guide_legend(override.aes = list(color = "black"))
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Genotype distribution by sample (sequence length)", x = NULL, y = "Sequence length (bp)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

# Plot 2: Violin plot of coverage distribution per sample
p2 <- ggplot(coverage %>% filter(coverage > 0), aes(x = sample, y = coverage, fill = sample)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "Coverage distribution by sample", 
    subtitle = "Violin plots with quartiles (log10 scale, zero coverage excluded)",
    x = NULL, 
    y = "Coverage (log10)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Plot 3: Genotypes across nodes (heatmap with rectangles scaled by node length)
p3 <- ggplot(zyg_long %>% filter(!is.na(genotype))) +
  geom_rect(
    aes(xmin = cum_start, xmax = cum_end, ymin = 0, ymax = 1, fill = genotype),
    color = NA
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(
    values = c("0" = "#2166ac", "0/0" = "#2166ac", "0/1" = "#fddbc7", 
               "1" = "#b2182b", "1/1" = "#b2182b",
               "." = "grey80", "./." = "grey80"),
    na.value = "grey90",
    name = "Genotype",
    guide = guide_legend(override.aes = list(color = "black"))
  ) +
  facet_wrap(~sample, ncol = 1, strip.position = "right") +
  labs(
    title = "Genotypes across nodes", 
    subtitle = "Rectangle width proportional to node length",
    x = "Panenomic position (bp)", 
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Plot 3_het: Heterozygous genotypes across nodes (only 0/1 cases)
# Check if there are any heterozygous genotypes
has_het <- any(zyg_long$genotype == "0/1", na.rm = TRUE)

if (has_het) {
  p3_het <- ggplot(zyg_long %>% filter(genotype == "0/1")) +
    geom_rect(
      aes(xmin = cum_start, xmax = cum_end, ymin = 0, ymax = 1),
      fill = "black",
    ) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    facet_wrap(~sample, ncol = 1, strip.position = "right") +
    labs(
      title = "Heterozygous regions (0/1) across nodes", 
      subtitle = "Highlighting only heterozygous genotypes with black borders; rectangle width proportional to node length",
      x = "Panenomic position (bp)", 
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.text.y = element_text(angle = 0),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Print heterozygous statistics
  het_stats <- zyg_long %>%
    filter(genotype == "0/1") %>%
    group_by(sample) %>%
    summarise(
      het_nodes = n(),
      het_length = sum(length, na.rm = TRUE),
      .groups = "drop"
    )
  
} else {
  # Create empty plot if no heterozygous genotypes
  p3_het <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No heterozygous (0/1) genotypes found", 
             size = 6, hjust = 0.5, vjust = 0.5) +
    labs(title = "Heterozygous regions (0/1) across nodes") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
}

# Plot 3b: Panplexity mask across nodes (only show low-complexity nodes)
# Check which categories exist in the data
panplexity_categories <- panplexity_mask %>% 
  filter(!is.na(masked)) %>% 
  pull(masked_label) %>% 
  unique()

# Build the plot
p3b <- ggplot(panplexity_mask %>% filter(!is.na(masked))) +
  geom_rect(
    aes(xmin = cum_start, xmax = cum_end, ymin = 0, ymax = 1, fill = masked_label),
    color = NA
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(
    values = c("Low-complexity" = "#d73027", "Not masked" = "transparent"),
    name = "Panplexity mask",
    breaks = c("Low-complexity", "Not masked"),
    drop = FALSE
  ) +
  facet_wrap(~"Pangenome", ncol = 1, strip.position = "right") +
  labs(
    title = "Panplexity mask across nodes", 
    subtitle = "Rectangle width proportional to node length (red = low-complexity nodes only)",
    x = "Panenomic position (bp)", 
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Override legend aesthetics only if needed (to show transparent as white)
if (length(panplexity_categories) > 1) {
  p3b <- p3b + guides(fill = guide_legend(
    override.aes = list(
      fill = c("Low-complexity" = "#d73027", "Not masked" = "white"),
      color = c("Low-complexity" = "black", "Not masked" = "black")
    )
  ))
}

# Plot 3c: Node coverage mask across nodes (only show masked nodes)
# Check which categories exist in the data
coverage_mask_categories <- node_coverage_mask %>% 
  filter(!is.na(coverage_masked)) %>% 
  pull(coverage_masked_label) %>% 
  unique()

# Build the plot
p3c <- ggplot(node_coverage_mask %>% filter(!is.na(coverage_masked))) +
  geom_rect(
    aes(xmin = cum_start, xmax = cum_end, ymin = 0, ymax = 1, fill = coverage_masked_label),
    color = NA
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(
    values = c("Masked" = "#e66101", "Not masked" = "transparent"),
    name = "Node coverage mask",
    breaks = c("Masked", "Not masked"),
    drop = FALSE
  ) +
  facet_wrap(~"Node coverage", ncol = 1, strip.position = "right") +
  labs(
    title = "Node coverage mask across nodes", 
    subtitle = "Rectangle width proportional to node length (orange = masked nodes only)",
    x = "Panenomic position (bp)", 
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Override legend aesthetics only if needed (to show transparent as white)
if (length(coverage_mask_categories) > 1) {
  p3c <- p3c + guides(fill = guide_legend(
    override.aes = list(
      fill = c("Masked" = "#e66101", "Not masked" = "white"),
      color = c("Masked" = "black", "Not masked" = "black")
    )
  ))
}

# Plot 4: Coverage across nodes (rectangles scaled by node length)
# Calculate max coverage per sample and where it occurs
max_cov_per_sample <- coverage %>%
  filter(!is.na(coverage), coverage > 0) %>%
  group_by(sample) %>%
  filter(coverage == max(coverage)) %>%
  slice(1) %>%  # In case there are multiple nodes with the same max coverage, take the first
  mutate(
    max_coverage = coverage,
    max_position = (cum_start + cum_end) / 2  # Use the center of the node
  ) %>%
  select(sample, max_coverage, max_position, node, cum_start, cum_end) %>%
  ungroup()

p4 <- ggplot(coverage %>% filter(!is.na(coverage), coverage > 0)) +
  geom_rect(
    aes(xmin = cum_start, xmax = cum_end, ymin = 0, ymax = coverage, fill = sample),
    color = NA, alpha = 0.8
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  # Add vertical lines at max coverage positions
  geom_vline(data = max_cov_per_sample, aes(xintercept = max_position), 
             color = "black", linetype = "dashed", linewidth = 0.6, alpha = 0.6) +
  # Add horizontal lines for max coverage values
  geom_hline(data = max_cov_per_sample, aes(yintercept = max_coverage), 
             color = "black", linetype = "dashed", linewidth = 0.6, alpha = 0.6) +
  # Add text annotation for max coverage value
  geom_text(data = max_cov_per_sample, 
            aes(x = Inf, y = max_coverage, label = sprintf("Max: %.1f", max_coverage)),
            hjust = 1.2, vjust = 1.5, size = 3.5, color = "black", fontface = "bold") +
  facet_wrap(~sample, ncol = 1, scales = "fixed", strip.position = "right") +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "Coverage across nodes", 
    subtitle = "Rectangle width proportional to node length (log10 scale, zero coverage excluded). Dashed lines show maximum coverage.",
    x = "Panenomic position (bp)", 
    y = "Coverage (log10)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text.y = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Plot 5: Low coverage across nodes (capped at median/2 for each sample)
# Calculate median coverage per sample for capping
median_cov_per_sample <- coverage %>%
  filter(!is.na(coverage), coverage > 0) %>%
  group_by(sample) %>%
  summarise(
    median_coverage = median(coverage),
    cap_coverage = median(coverage) / 5,
    .groups = "drop"
  )

# Prepare capped coverage data
coverage_capped <- coverage %>%
  left_join(median_cov_per_sample, by = "sample") %>%
  mutate(coverage_capped = pmin(coverage, cap_coverage, na.rm = TRUE))

p5 <- ggplot(coverage_capped %>% filter(!is.na(coverage_capped), coverage_capped > 0)) +
  geom_rect(
    aes(xmin = cum_start, xmax = cum_end, ymin = 0, ymax = coverage_capped, fill = sample),
    color = NA, alpha = 0.8
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  facet_wrap(~sample, ncol = 1, scales = "free_y", strip.position = "right") +
  labs(
    title = "Low coverage regions across nodes", 
    subtitle = "Y-axis capped at median/5 to highlight low coverage regions. Zero coverage excluded.",
    x = "Panenomic position (bp)", 
    y = "Coverage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text.y = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Additional plot: Coverage vs Genotype
p6 <- ggplot(
  zyg_long %>% 
    left_join(coverage %>% select(sample, node, coverage), by = c("sample", "node")) %>%
    filter(!is.na(gt_num), !is.na(coverage))
) +
  geom_violin(aes(x = factor(gt_num), y = coverage, fill = factor(gt_num)), alpha = 0.7) +
  geom_boxplot(aes(x = factor(gt_num), y = coverage), width = 0.1, outlier.shape = NA) +
  scale_y_continuous(trans = "log10") +
  scale_x_discrete(labels = c("0" = "0/0", "1" = "0/1", "2" = "1/1")) +
  scale_fill_manual(
    values = c("0" = "#2166ac", "1" = "#fddbc7", "2" = "#b2182b"),
    labels = c("0/0", "0/1", "1/1"),
    name = "Genotype"
  ) +
  facet_wrap(~sample, scales = "free_y") +
  labs(
    title = "Coverage distribution by genotype call",
    subtitle = "Violin plots showing coverage for each genotype (log10 scale)",
    x = "Genotype", 
    y = "Coverage (log10)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save plots
cat("Saving plots to:", output_pdf, "\n")
pdf(output_pdf, width = 16, height = 36)

# Layout with all plots including the new heterozygous plot
print(
  (p1 | p1b) / 
    (p2 | p6) / 
    p3 / 
    p3_het /
    p3b /
    p3c /
    p4 /
    p5 + 
    plot_layout(heights = c(2, 2, 1.5, 1.5, 0.2, 0.2, 2.5, 2.5))
)

dev.off()

cat("\n=== Summary Statistics ===\n")
cat("\nCoverage Statistics:\n")
print(cov_stats, n = Inf)

cat("\nZero Coverage Nodes:\n")
zero_cov <- coverage %>%
  group_by(sample) %>%
  summarise(
    zero_nodes = sum(coverage == 0),
    total_nodes = n(),
    pct_zero = (zero_nodes / total_nodes) * 100,
    .groups = "drop"
  )
print(zero_cov, n = Inf)

cat("\nGenotype Counts:\n")
genotype_summary <- zyg_stats %>% 
  pivot_wider(names_from = genotype_label, values_from = count, values_fill = 0)
print(genotype_summary, n = Inf)

cat("\nGenotype Percentages:\n")
genotype_pct <- zyg_stats_pct %>% 
  select(sample, genotype_label, percentage) %>%
  pivot_wider(names_from = genotype_label, values_from = percentage, values_fill = 0)
print(genotype_pct, n = Inf)

# Print sequence length statistics
cat("\nGenotype Sequence Lengths (bp):\n")
genotype_len_summary <- zyg_stats_len %>% 
  pivot_wider(names_from = genotype_label, values_from = total_length, values_fill = 0)
print(genotype_len_summary, n = Inf)

cat("\nGenotype Sequence Length Percentages:\n")
genotype_len_pct <- zyg_stats_len_pct %>% 
  select(sample, genotype_label, percentage) %>%
  pivot_wider(names_from = genotype_label, values_from = percentage, values_fill = 0)
print(genotype_len_pct, n = Inf)

# Print heterozygous statistics if they exist
if (has_het) {
  cat("\nHeterozygous (0/1) Statistics:\n")
  print(het_stats, n = Inf)
  
  # Calculate percentage of heterozygous nodes
  het_pct <- zyg_long %>%
    group_by(sample) %>%
    summarise(
      total_nodes = n(),
      het_nodes = sum(genotype == "0/1", na.rm = TRUE),
      het_percentage = (het_nodes / total_nodes) * 100,
      .groups = "drop"
    )
  
  cat("\nHeterozygous Node Percentages:\n")
  print(het_pct, n = Inf)
  
  # Calculate percentage of heterozygous sequence length
  het_len_pct <- zyg_long %>%
    group_by(sample) %>%
    summarise(
      total_length = sum(length, na.rm = TRUE),
      het_length = sum(length[genotype == "0/1"], na.rm = TRUE),
      het_length_percentage = (het_length / total_length) * 100,
      .groups = "drop"
    )
  
  cat("\nHeterozygous Sequence Length Percentages:\n")
  print(het_len_pct, n = Inf)
} else {
  cat("\nNo heterozygous (0/1) genotypes found in the data.\n")
}

# Print panplexity mask statistics
cat("\nPanplexity Mask Statistics:\n")
mask_stats <- panplexity_mask %>%
  group_by(masked_label) %>%
  summarise(
    node_count = n(),
    total_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    node_pct = (node_count / sum(node_count)) * 100,
    length_pct = (total_length / sum(total_length)) * 100
  )
print(mask_stats, n = Inf)

# Print node coverage mask statistics
cat("\nNode Coverage Mask Statistics:\n")
coverage_mask_stats <- node_coverage_mask %>%
  group_by(coverage_masked_label) %>%
  summarise(
    node_count = n(),
    total_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    node_pct = (node_count / sum(node_count)) * 100,
    length_pct = (total_length / sum(total_length)) * 100
  )
print(coverage_mask_stats, n = Inf)

# Print maximum coverage node information
cat("\nMaximum Coverage Node Information:\n")
max_cov_info <- max_cov_per_sample %>%
  select(sample, node, max_coverage, max_position) %>%
  mutate(
    position_desc = sprintf("Node %d at position %.0f bp", node, max_position)
  )
print(max_cov_info, n = Inf)

cat("\nPlots saved to:", output_pdf, "\n")

