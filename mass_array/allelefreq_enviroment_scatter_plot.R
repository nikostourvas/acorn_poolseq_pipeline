library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(forcats)
library(viridis)

# Read the SNP positions file
snp_positions <- fread(
  "2025-08-25_massarray_proposed_snps_bio1-3-10-12-15-18_q0001-bf10_variable-minMAF.tsv"
)

# Clean chromosome column and create SNP ID
snp_positions_clean <- snp_positions |>
  mutate(
    chr = as.integer(gsub("chr", "", chr))
  ) |>
  mutate(
    chr = sprintf("%02d", chr),
    chrom_pos = paste0("Qrob_Chr", chr, "_", pos)
  )

print(paste("Loaded", nrow(snp_positions_clean), "SNP positions"))
print("Bioclimate variables in dataset:")
print(table(snp_positions_clean$bioclim))

# Load allele frequency data
# You'll need to replace this with your actual allele frequency file
# Expected format: columns for chr, pos, and population allele frequencies
afs <- fread(
  "/home/geneticsShare/LFMM_ANALYSES/dat/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_AlleleFrequencyTable.txt",
  header = TRUE
)

# Keep only the chrom_pos rows present in snp_positions_clean
dataset <- left_join(snp_positions_clean, afs)

# release afs from memory
rm(afs)

# Transform to tidy format:
dataset <- dataset |>
  pivot_longer(
    cols = -c(chr, pos, MAF, bioclim, method, chrom_pos),
    names_to = "population",
    values_to = "allele_freq"
  )

# Load environmental data
env_data <- fread(
  "/home/geneticsShare/LFMM_ANALYSES/dat/20240805_ACORN_dem_TOPO_by_POP_NT.csv"
)
# Only keep the columns corresponding to the bioclimate variables in snp_positions_clean (bioclim)
env_data$Plot_ID <- as.character(env_data$Plot_ID)
env_data <- env_data[, c(1, 18:36)]
env_data <- pivot_longer(
  env_data,
  cols = -Plot_ID,
  names_to = "bioclim",
  values_to = "env_value"
)

env_data <- env_data |>
  # select the rows with Plot_ID in population names
  filter(Plot_ID %in% dataset$population) |>
  filter(bioclim %in% unique(dataset$bioclim))

# Merge allele frequency data with environmental data
dataset <- full_join(
  dataset,
  env_data,
  by = c("population" = "Plot_ID", "bioclim" = "bioclim")
)

# Create scatter plots only for bioclimate variables present in the SNP data
create_targeted_scatter_plots <- function(
  data,
  output_dir = "/home/geneticsShare/plots/"
) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Get only the bioclimate variables that are actually present in the SNP data
  present_bioclim_vars <- unique(data$bioclim)

  message(paste(
    "Creating plots for",
    length(present_bioclim_vars),
    "bioclimate variables:"
  ))
  message(paste(present_bioclim_vars, collapse = ", "))

  for (bioclim_var in present_bioclim_vars) {
    # Filter data for current bioclimate variable
    bioclim_data <- data |>
      filter(bioclim == bioclim_var)

    if (nrow(bioclim_data) == 0) {
      message(paste("No data for", bioclim_var, "- skipping"))
      next
    }

    # Count unique SNPs for this bioclimate variable
    n_snps_total <- length(unique(bioclim_data$chrom_pos))
    message(paste("Processing", bioclim_var, "with", n_snps_total, "SNPs"))

    # Order SNPs by MAF for consistent plotting
    bioclim_data <- bioclim_data |>
      mutate(chrom_pos = fct_reorder(chrom_pos, MAF, .desc = TRUE))

    # Create the plot
    p <- ggplot(bioclim_data, aes(x = env_value, y = allele_freq)) +
      geom_point(alpha = 0.7, color = "steelblue", size = 1.5) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "red",
        alpha = 0.3,
        linewidth = 0.8
      ) +
      facet_wrap(~chrom_pos, scales = "free_x", ncol = 5) +
      scale_y_continuous(
        limits = c(0, 1),
        labels = scales::percent_format(accuracy = 1)
      ) +
      labs(
        title = paste("Allele Frequency vs", bioclim_var),
        subtitle = paste(
          "Showing",
          length(unique(bioclim_data$chrom_pos)),
          "SNPs"
        ),
        x = paste(bioclim_var, "Environmental Value"),
        y = "Allele Frequency"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 7),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9),
        panel.grid.minor = element_blank()
      )

    # Calculate plot dimensions based on number of SNPs
    n_snps <- length(unique(bioclim_data$chrom_pos))
    n_cols <- 5
    n_rows <- ceiling(n_snps / n_cols)

    plot_width <- min(30, max(12, n_cols * 2.5))
    plot_height <- max(8, n_rows * 2.5)

    # Save the plot
    filename <- file.path(output_dir, paste0(bioclim_var, "_scatter_plots.pdf"))
    ggsave(
      filename,
      p,
      width = plot_width,
      height = plot_height,
      units = "in",
      limitsize = FALSE
    )

    message(paste("  Saved plot:", filename))
  }

  message(paste(
    "Completed scatter plots for",
    length(present_bioclim_vars),
    "bioclimate variables"
  ))
}

# Function to create a summary plot showing all bioclimate variables
create_summary_plot <- function(
  data,
  output_dir = "/home/geneticsShare/plots/"
) {
  # Create a summary of correlations by bioclimate variable
  correlation_summary <- data |>
    group_by(bioclim, chrom_pos) |>
    summarise(
      correlation = cor(env_value, allele_freq, use = "complete.obs"),
      n_points = n(),
      .groups = "drop"
    ) |>
    filter(!is.na(correlation))

  # Plot distribution of correlations by bioclimate variable
  p_summary <- ggplot(
    correlation_summary,
    aes(x = bioclim, y = correlation, fill = bioclim)
  ) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    scale_fill_viridis_d() +
    labs(
      title = "Distribution of Allele Frequency-Environment Correlations",
      subtitle = "Each point represents one SNP",
      x = "Bioclimate Variable",
      y = "Correlation (Allele Freq ~ Env Value)",
      fill = "Bioclimate"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  filename_summary <- file.path(output_dir, "correlation_summary.pdf")
  ggsave(filename_summary, p_summary, width = 10, height = 6, units = "in")

  message(paste("Saved summary plot:", filename_summary))
}

# Example usage (uncomment when you have the data):
# create_targeted_scatter_plots(plot_data)
# create_summary_plot(plot_data)

# Print summary of what will be plotted
print("\nSNPs per bioclimate variable:")
snp_summary <- snp_positions_clean |>
  count(bioclim, name = "n_snps") |>
  arrange(desc(n_snps))
print(snp_summary)

print(paste(
  "\nThis script will create scatter plots for",
  nrow(snp_summary),
  "bioclimate variables"
))
print("Each plot will show allele frequency vs environmental values")
print("All SNPs will be included in the plots (no limit applied)")

print("\nTo complete the analysis, you need to:")
print("1. Load your allele frequency data file")
print("2. Load your environmental data file")
print("3. Merge the datasets")
print("4. Call create_targeted_scatter_plots() with the merged data")
