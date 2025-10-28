#!/usr/bin/env Rscript

# Script to convert outlier CSV files to BED format
#
# This script processes CSV files containing outlier genomic windows
# and converts them to BED format files.
#
# It calculates window coordinates based on the 'index' column (e.g., "1_window_00976")
# and uses the 'envfactor' column as the name field in the BED file.
#
# Usage:
#   Rscript wza_res_to_bed.R [base_dir] [pattern] [window_size]
#
# Examples:
#   # Process all .csv files in current directory
#   Rscript wza_res_to_bed.R
#
#   # Process files in specific directory
#   Rscript wza_res_to_bed.R results/my_csvs
#
#   # Custom pattern and window size
#   Rscript wza_res_to_bed.R . "my_data_.*\\.csv$" 10000
#
# Window name format (in 'index' column): chr_window_number (e.g., "2_window_01619")
#
# Output: BED files with columns: chrom, chromStart, chromEnd, name (from envfactor)
#
# Author: Adapted from a script by Nikolaos Tourvas
# Date: 2025-10-28

# Set options to avoid scientific notation in output
options(scipen = 500)

library(dplyr)
library(readr)

#' Convert window name to genomic coordinates
#' (Reused from original script)
#' @param window_name Character string in format "chr_window_number"
#' @param window_size Size of each window in base pairs (default: 5000)
#' @return Data frame with chrom, chromStart, chromEnd, window_name
parse_window_to_bed <- function(window_name, window_size = 5000) {
  # Split the window name by underscore
  parts <- strsplit(window_name, "_")[[1]]
  
  if (length(parts) != 3 || parts[2] != "window") {
    warning(paste("Invalid window name format:", window_name))
    return(NULL)
  }
  
  # Extract chromosome and window number
  chrom <- paste0("chr", parts[1])
  window_num <- as.numeric(parts[3])
  
  # Calculate 0-based coordinates (BED format: start inclusive, end exclusive)
  # Window number is 1-based, so window 1 starts at position 0
  chromStart <- (window_num - 1) * window_size
  chromEnd <- chromStart + window_size
  
  return(data.frame(
    chrom = chrom,
    chromStart = chromStart,
    chromEnd = chromEnd,
    window_name = window_name, # Original index, used for joining
    stringsAsFactors = FALSE
  ))
}

#' Process a single outlier CSV file and generate BED file
#' (Modified from original script to handle CSV input)
#' @param file_path Path to the input CSV file
#' @param window_size Size of each window in base pairs (default: 5000)
process_csv_file <- function(file_path, window_size = 5000) {
  cat("Processing:", file_path, "\n")
  
  # Check if file exists
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(FALSE)
  }
  
  tryCatch({
    # Read the CSV file
    csv_data <- read_csv(file_path, col_types = cols())
    
    if (nrow(csv_data) == 0) {
      warning(paste("CSV file is empty:", file_path))
      return(FALSE)
    }
    
    # Check for required columns
    required_cols <- c("index", "envfactor")
    if (!all(required_cols %in% names(csv_data))) {
      missing_cols <- required_cols[!required_cols %in% names(csv_data)]
      warning(paste("Missing required columns:", paste(missing_cols, collapse=", "), "in", file_path))
      return(FALSE)
    }

    # Convert each window name in the 'index' column to BED format
    bed_list <- lapply(csv_data$index, parse_window_to_bed, window_size = window_size)
    
    # Identify and remove NULL entries (invalid window names)
    valid_indices <- !sapply(bed_list, is.null)
    bed_list <- bed_list[valid_indices]
    
    if (length(bed_list) == 0) {
      warning(paste("No valid windows could be parsed from:", file_path))
      return(FALSE)
    }
    
    # Combine into a single data frame
    bed_df <- do.call(rbind, bed_list)
    
    # Add the 'envfactor' as the name column
    # Filter original csv_data to match only validly parsed windows
    bed_df$name <- csv_data$envfactor[valid_indices]
    
    # Select and reorder columns for final BED format
    bed_df <- bed_df %>%
      select(chrom, chromStart, chromEnd, name)
    
    # Sort by chromosome and start position
    bed_df <- bed_df[order(bed_df$chrom, bed_df$chromStart), ]
    
    # (Optional: kept from original script)
    bed_df$chrom <- sub("^chrQrob\\.", "", bed_df$chrom)
    
    # Ensure coordinates are formatted as integers (no scientific notation)
    bed_df$chromStart <- as.integer(bed_df$chromStart)
    bed_df$chromEnd <- as.integer(bed_df$chromEnd)
    
    # Generate output file name (from .csv to .bed)
    output_file <- sub("\\.csv$", ".bed", file_path, ignore.case = TRUE)
    
    # Write BED file using cat to ensure proper formatting
    bed_lines <- paste(bed_df$chrom, bed_df$chromStart, bed_df$chromEnd, bed_df$name, sep = "\t")
    writeLines(bed_lines, output_file)
    
    cat("Generated BED file:", output_file, "with", nrow(bed_df), "regions\n")
    return(TRUE)
    
  }, error = function(e) {
    warning(paste("Error processing", file_path, ":", e$message))
    return(FALSE)
  })
}

#' Find and process all outlier files
#' (Modified to change default pattern and call new function)
#' @param base_dir Base directory to search for files (default: current directory)
#' @param pattern Pattern to match file names (default: ".*\\.csv$")
#' @param window_size Size of each window in base pairs (default: 5000)
main <- function(base_dir = ".", pattern = ".*\\.csv$", window_size = 5000) {
  cat("Searching for outlier CSV files in:", base_dir, "\n")
  cat("Pattern:", pattern, "\n")
  cat("Window size:", window_size, "bp\n\n")
  
  # Find all matching files recursively
  all_files <- list.files(base_dir, 
                         pattern = pattern, 
                         recursive = TRUE, 
                         full.names = TRUE)
  
  if (length(all_files) == 0) {
    cat("No files found matching pattern:", pattern, "\n")
    return()
  }
  
  cat("Found", length(all_files), "files to process:\n")
  for (file in all_files) {
    cat(" -", file, "\n")
  }
  cat("\n")
  
  # Process each file
  success_count <- 0
  for (file_path in all_files) {
    # Call the new CSV processing function
    if (process_csv_file(file_path, window_size)) {
      success_count <- success_count + 1
    }
  }
  
  cat("\nSummary:\n")
  cat("Total files found:", length(all_files), "\n")
  cat("Successfully processed:", success_count, "\n")
  cat("Failed:", length(all_files) - success_count, "\n")
}

# Command line argument parsing
# (Reused from original, updated help text)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # Default: search in current directory
  main()
} else if (length(args) == 1) {
  # Custom base directory
  main(base_dir = args[1])
} else if (length(args) == 2) {
  # Custom base directory and pattern
  main(base_dir = args[1], pattern = args[2])
} else if (length(args) == 3) {
  # Custom base directory, pattern, and window size
  main(base_dir = args[1], pattern = args[2], window_size = as.numeric(args[3]))
} else {
  cat("Usage: Rscript convert_csv_to_bed.R [base_dir] [pattern] [window_size]\n")
  cat("  base_dir: Directory to search for files (default: current directory)\n")
  cat("  pattern: File pattern to match (default: '.*\\\\.csv$')\n")
  cat("  window_size: Size of genomic windows in bp (default: 5000)\n")
  quit(status = 1)
}