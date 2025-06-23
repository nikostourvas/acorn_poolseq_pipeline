# ---- Three-Way Table Comparison Function ----
# Purpose: Compare three tables based on "index" column and generate overlap analysis with Venn diagrams

# ---- Load required libraries ----
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(data.table)
library(viridis)
library(VennDiagram)
library(grid)
library(gridExtra)

# ---- Define color palette ----
my_pal <- c(
  neutral_light = "#BBBBBB", neutral_dark = "#000000",
  blue_light    = "#6699CC", blue_dark    = "#004488",
  gold_light    = "#EECC66", gold_dark    = "#997700",
  sig_red       = "#CC3311"
)

# ---- Venn diagram function (re-adjusted) ----
generate_venn <- function(set, name, output_dir = ".") {
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Full path for the output file
  full_path <- file.path(output_dir, paste0(name, ".svg"))
  
  svg(full_path, width = 6, height = 6, family = "sans")
  
  venn.grob <- venn.diagram(
    x              = set,
    filename       = NULL,          # keep NULL – we draw it ourselves
    lwd            = 2,
    fill           = c(my_pal[[1]], my_pal[[3]], my_pal[[5]]),
    ## region-size numbers
    cex            = 1.2,
    fontface       = "plain",
    fontfamily     = "sans",        # <=  region labels
    print.mode     = "raw",         # "raw", "percent"
    ## category (set) labels
    cat.cex        = 1.0,
    cat.fontface   = "bold",
    cat.fontfamily = "sans",        # <=  set labels
    cat.dist       = 0.06,
    cat.pos        = 0,
    margin         = 0.04
  )
  
  grid.draw(venn.grob)
  dev.off()        # closes the SVG file
  
  message("✓ Venn diagram saved to ", full_path)
}

# ---- Helper function to read table with index column ----
read_table_with_index <- function(file_path, table_name = NULL) {
  # Try different file formats
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    data <- fread(file_path)
  } else if (grepl("\\.txt$", file_path, ignore.case = TRUE)) {
    data <- fread(file_path)
  } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    data <- fread(file_path, sep = "\t")
  } else {
    # Try to read as CSV by default
    data <- fread(file_path)
  }
  
  # Check if "index" column exists
  if (!"index" %in% colnames(data)) {
    stop("Error: Column 'index' not found in file: ", file_path)
  }
  
  # Add table name if provided
  if (!is.null(table_name)) {
    data$table_source <- table_name
  }
  
  return(data)
}

# ---- Main three-way comparison function ----
compare_three_tables <- function(table1_path, table2_path, table3_path, 
                                table1_name = "Table1", table2_name = "Table2", table3_name = "Table3",
                                output_dir = "three_way_comparison_results",
                                generate_venn_diagram = TRUE,
                                n_windows = NULL) {
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Check if n_windows is provided for hypergeometric test
  if (is.null(n_windows)) {
    stop("Error: n_windows (total universe size) must be provided as an argument")
  }
  
  message("Starting three-way table comparison analysis...")
  message("Table 1: ", table1_name, " (", table1_path, ")")
  message("Table 2: ", table2_name, " (", table2_path, ")")
  message("Table 3: ", table3_name, " (", table3_path, ")")
  message("Total universe size (n_windows): ", n_windows)
  
  # Read the three tables
  message("\nReading tables...")
  
  tryCatch({
    table1 <- read_table_with_index(table1_path, table1_name)
    message("✓ ", table1_name, " loaded: ", nrow(table1), " rows")
  }, error = function(e) {
    stop("Error reading ", table1_name, ": ", e$message)
  })
  
  tryCatch({
    table2 <- read_table_with_index(table2_path, table2_name)
    message("✓ ", table2_name, " loaded: ", nrow(table2), " rows")
  }, error = function(e) {
    stop("Error reading ", table2_name, ": ", e$message)
  })
  
  tryCatch({
    table3 <- read_table_with_index(table3_path, table3_name)
    message("✓ ", table3_name, " loaded: ", nrow(table3), " rows")
  }, error = function(e) {
    stop("Error reading ", table3_name, ": ", e$message)
  })
  
  # Extract unique indices from each table
  message("\nExtracting unique indices...")
  indices1 <- unique(table1$index)
  indices2 <- unique(table2$index)
  indices3 <- unique(table3$index)
  
  message("✓ ", table1_name, ": ", length(indices1), " unique indices")
  message("✓ ", table2_name, ": ", length(indices2), " unique indices")
  message("✓ ", table3_name, ": ", length(indices3), " unique indices")
  
  # Calculate all possible overlaps
  message("\nCalculating overlaps...")
  
  # Pairwise overlaps
  overlap_1_2 <- intersect(indices1, indices2)
  overlap_1_3 <- intersect(indices1, indices3)
  overlap_2_3 <- intersect(indices2, indices3)
  
  # Three-way overlap
  overlap_1_2_3 <- intersect(intersect(indices1, indices2), indices3)
  
  # Unique to each set (not in any other)
  unique_1 <- setdiff(setdiff(indices1, indices2), indices3)
  unique_2 <- setdiff(setdiff(indices2, indices1), indices3)
  unique_3 <- setdiff(setdiff(indices3, indices1), indices2)
  
  # Create comprehensive overlap summary
  overlap_summary <- tibble(
    comparison = c(
      paste0(table1_name, " only"),
      paste0(table2_name, " only"),
      paste0(table3_name, " only"),
      paste0(table1_name, " ∩ ", table2_name, " only"),
      paste0(table1_name, " ∩ ", table3_name, " only"),
      paste0(table2_name, " ∩ ", table3_name, " only"),
      paste0(table1_name, " ∩ ", table2_name, " ∩ ", table3_name)
    ),
    count = c(
      length(unique_1),
      length(unique_2),
      length(unique_3),
      length(setdiff(overlap_1_2, overlap_1_2_3)),
      length(setdiff(overlap_1_3, overlap_1_2_3)),
      length(setdiff(overlap_2_3, overlap_1_2_3)),
      length(overlap_1_2_3)
    ),
    percentage = round(c(
      length(unique_1) / length(indices1) * 100,
      length(unique_2) / length(indices2) * 100,
      length(unique_3) / length(indices3) * 100,
      length(setdiff(overlap_1_2, overlap_1_2_3)) / length(indices1) * 100,
      length(setdiff(overlap_1_3, overlap_1_2_3)) / length(indices1) * 100,
      length(setdiff(overlap_2_3, overlap_1_2_3)) / length(indices2) * 100,
      length(overlap_1_2_3) / length(indices1) * 100
    ), 2)
  )
  
  # Create pairwise comparison table (similar to original script)
  pairwise_comparisons <- tibble(
    set_A = c(table1_name, table1_name, table2_name),
    set_B = c(table2_name, table3_name, table3_name),
    n_A = c(length(indices1), length(indices1), length(indices2)),
    n_B = c(length(indices2), length(indices3), length(indices3)),
    n_AB = c(length(overlap_1_2), length(overlap_1_3), length(overlap_2_3)),
    prop_AB = c(
      length(overlap_1_2) / min(length(indices1), length(indices2)),
      length(overlap_1_3) / min(length(indices1), length(indices3)),
      length(overlap_2_3) / min(length(indices2), length(indices3))
    )
  )
  
  # Calculate overlap of table 1 and table 2, then find overlap with table 3
  message("Calculating overlap of ", table1_name, " and ", table2_name, ", then overlap with ", table3_name)
  overlap_1_2 <- intersect(indices1, indices2)
  overlap_intersection_1_2_with_3 <- intersect(overlap_1_2, indices3)
  
  # Create special comparison for intersection vs table 3 with hypergeometric p-value
  # Following the pattern from pairwise_overlap_analysis.R
  union_comparison <- tibble(
    comparison_type = paste0("Overlap(", table1_name, ", ", table2_name, ") vs ", table3_name),
    n_overlap_1_2 = length(overlap_1_2),
    n_table3 = length(indices3),
    n_overlap = length(overlap_intersection_1_2_with_3),
    prop_overlap_of_intersection = length(overlap_intersection_1_2_with_3) / length(overlap_1_2),
    prop_overlap_of_table3 = length(overlap_intersection_1_2_with_3) / length(indices3),
    prop_AB = length(overlap_intersection_1_2_with_3) / min(length(overlap_1_2), length(indices3)),
    pvalue = phyper(q = length(overlap_intersection_1_2_with_3) - 1,
                   m = length(overlap_1_2),
                   n = n_windows - length(overlap_1_2),
                   k = length(indices3),
                   lower.tail = FALSE)
  )
  
  # Save results
  message("\nSaving results...")
  
  # Save overlap summary
  overlap_summary_file <- file.path(output_dir, "overlap_summary.csv")
  write_csv(overlap_summary, overlap_summary_file)
  message("✓ Overlap summary saved to ", overlap_summary_file)
  
  # Save pairwise comparisons
  pairwise_file <- file.path(output_dir, "pairwise_comparisons.csv")
  write_csv(pairwise_comparisons, pairwise_file)
  message("✓ Pairwise comparisons saved to ", pairwise_file)
  
  # Save union comparison
  intersection_comparison_file <- file.path(output_dir, "intersection_comparison.csv")
  write_csv(union_comparison, intersection_comparison_file)
  message("✓ Intersection comparison saved to ", intersection_comparison_file)
  
  # Save individual overlap sets - only three-way overlap
  write_lines(overlap_1_2_3, file.path(output_dir, paste0("overlap_all_three_", table1_name, "_", table2_name, "_", table3_name, ".txt")))
  
  message("✓ Three-way overlap file saved")
  
  # Generate Venn diagram if requested
  if (generate_venn_diagram) {
    message("\nGenerating Venn diagram...")
    
    # Prepare data for Venn diagram
    venn_data <- list(
      indices1,
      indices2,
      indices3
    )
    names(venn_data) <- c(table1_name, table2_name, table3_name)
    
    # Generate Venn diagram
    venn_filename <- paste0("venn_diagram_", table1_name, "_", table2_name, "_", table3_name)
    generate_venn(venn_data, venn_filename, output_dir)
  }
  
  # Create summary visualization
  message("\nCreating summary visualizations...")
  
  # Overlap counts bar plot
  overlap_plot <- ggplot(overlap_summary, aes(x = reorder(comparison, count), y = count)) +
    geom_col(fill = my_pal[["blue_dark"]], alpha = 0.8) +
    geom_text(aes(label = count), hjust = -0.1, size = 3) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Overlap Analysis: Number of Shared Indices",
      subtitle = paste("Comparison of", table1_name, ",", table2_name, "and", table3_name),
      x = "Comparison Type",
      y = "Number of Indices"
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave(file.path(output_dir, "overlap_counts_barplot.png"), 
         overlap_plot, width = 10, height = 6, dpi = 300)
  
  # Pairwise heatmap
  pairwise_plot <- ggplot(pairwise_comparisons, aes(x = set_A, y = set_B, fill = prop_AB)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = paste0(round(prop_AB * 100, 1), "%")), color = "white", size = 4) +
    scale_fill_viridis_c(direction = -1, name = "Proportion\nOverlap") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Pairwise Overlap Proportions",
      subtitle = "Proportion of smaller set that overlaps with larger set",
      x = "Set A",
      y = "Set B"
    )
  
  ggsave(file.path(output_dir, "pairwise_overlap_heatmap.png"), 
         pairwise_plot, width = 8, height = 6, dpi = 300)
  
  message("✓ Summary visualizations saved")
  
  # Print summary to console
  message("\n==== SUMMARY RESULTS ====")
  print(overlap_summary)
  message("\n==== PAIRWISE COMPARISONS ====")
  print(pairwise_comparisons)
  message("\n==== INTERSECTION COMPARISON ====")
  print(union_comparison)
  
  # Return results as a list
  results <- list(
    overlap_summary = overlap_summary,
    pairwise_comparisons = pairwise_comparisons,
    union_comparison = union_comparison,
    indices = list(
      table1 = indices1,
      table2 = indices2,
      table3 = indices3
    ),
    overlaps = list(
      table1_table2 = overlap_1_2,
      table1_table3 = overlap_1_3,
      table2_table3 = overlap_2_3,
      all_three = overlap_1_2_3,
      intersection_1_2_with_3 = overlap_intersection_1_2_with_3
    ),
    unique_sets = list(
      table1_only = unique_1,
      table2_only = unique_2,
      table3_only = unique_3
    ),
    intersection_1_2 = overlap_1_2
  )
  
  message("\n✓ Analysis complete! Results saved to: ", output_dir)
  return(results)
}

# ---- Example usage function ----
# run_example_comparison <- function() {
#   # Example usage - replace with your actual file paths
#   example_results <- compare_three_tables(
#     table1_path = "path/to/your/first_table.csv",
#     table2_path = "path/to/your/second_table.csv", 
#     table3_path = "path/to/your/third_table.csv",
#     table1_name = "Dataset_A",
#     table2_name = "Dataset_B", 
#     table3_name = "Dataset_C",
#     output_dir = "example_comparison_results",
#     generate_venn_diagram = TRUE,
#     n_windows = 122034  # Total number of windows/features in the universe
#   )
  
#   return(example_results)
# }

# ---- Print usage instructions ----
message("Three-way table comparison function loaded!")
message("Usage:")
message("  results <- compare_three_tables(")
message("    table1_path = 'path/to/table1.csv',")
message("    table2_path = 'path/to/table2.csv',")
message("    table3_path = 'path/to/table3.csv',")
message("    table1_name = 'Table1_Name',")
message("    table2_name = 'Table2_Name',")
message("    table3_name = 'Table3_Name',")
message("    output_dir = 'output_directory',")
message("    generate_venn_diagram = TRUE,")
message("    n_windows = 122034  # Required: total number of windows/features")
message("  )")
