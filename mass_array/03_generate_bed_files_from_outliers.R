#!/usr/bin/env Rscript

# Script to convert outlier genomic windows to BED format
# 
# This script processes text files containing outlier genomic window names
# and converts them to BED format files with genomic coordinates.
#
# Usage:
#   Rscript generate_bed_files_from_outliers.R [base_dir] [pattern] [window_size]
#
# Examples:
#   # Process all overlap_*.txt files in current directory
#   Rscript generate_bed_files_from_outliers.R
#
#   # Process files in specific directory
#   Rscript generate_bed_files_from_outliers.R results/comparisons_wp3_wp4/Paired_GEA
#
#   # Custom pattern and window size
#   Rscript generate_bed_files_from_outliers.R . "outlier_.*\\.txt$" 10000
#
# Window name format: chr_window_number (e.g., "2_window_01619")
# - chr: chromosome number or scaffold name
# - window_number: 1-based window number (zero-padded to 5 digits)
# - Window size: default 5000 bp (5kb windows)
#
# Output: BED files with columns: chrom, chromStart, chromEnd, window_name
# 
# Author: Nikolaos Tourvas
# Date: 2025-06-20

# Set options to avoid scientific notation in output
# https://divingintogeneticsandgenomics.com/post/three-gotchas-when-using-r-for-genomic-data-analysis/
options(scipen = 500)

library(dplyr)
library(readr)

#' Convert window name to genomic coordinates
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
    window_name = window_name,
    stringsAsFactors = FALSE
  ))
}

#' Process a single outlier file and generate BED file
#' @param file_path Path to the input text file
#' @param window_size Size of each window in base pairs (default: 5000)
process_outlier_file <- function(file_path, window_size = 5000) {
  cat("Processing:", file_path, "\n")
  
  # Check if file exists
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(FALSE)
  }
  
  # Read the window names
  tryCatch({
    window_names <- readLines(file_path)
    # Remove empty lines and whitespace
    window_names <- trimws(window_names)
    window_names <- window_names[window_names != ""]
    
    if (length(window_names) == 0) {
      warning(paste("No valid window names found in:", file_path))
      return(FALSE)
    }
    
    # Convert each window name to BED format
    bed_list <- lapply(window_names, parse_window_to_bed, window_size = window_size)
    
    # Remove NULL entries (invalid window names)
    bed_list <- bed_list[!sapply(bed_list, is.null)]
    
    if (length(bed_list) == 0) {
      warning(paste("No valid windows could be parsed from:", file_path))
      return(FALSE)
    }
    
    # Combine into a single data frame
    bed_df <- do.call(rbind, bed_list)
    
    # Sort by chromosome and start position
    bed_df <- bed_df[order(bed_df$chrom, bed_df$chromStart), ]
    
    bed_df$chrom <- sub("^chrQrob\\.", "", bed_df$chrom)
    
    # Ensure coordinates are formatted as integers (no scientific notation)
    bed_df$chromStart <- as.integer(bed_df$chromStart)
    bed_df$chromEnd <- as.integer(bed_df$chromEnd)
    
    # Generate output file name
    output_file <- sub("\\.txt$", ".bed", file_path)
    
    # Write BED file using cat to ensure proper formatting
    bed_lines <- paste(bed_df$chrom, bed_df$chromStart, bed_df$chromEnd, bed_df$window_name, sep = "\t")
    writeLines(bed_lines, output_file)
    
    cat("Generated BED file:", output_file, "with", nrow(bed_df), "regions\n")
    return(TRUE)
    
  }, error = function(e) {
    warning(paste("Error processing", file_path, ":", e$message))
    return(FALSE)
  })
}

#' Find and process all outlier files
#' @param base_dir Base directory to search for files (default: current directory)
#' @param pattern Pattern to match file names (default: "overlap_.*\\.txt$")
#' @param window_size Size of each window in base pairs (default: 5000)
main <- function(base_dir = ".", pattern = "overlap_.*\\.txt$", window_size = 5000) {
  cat("Searching for outlier files in:", base_dir, "\n")
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
    if (process_outlier_file(file_path, window_size)) {
      success_count <- success_count + 1
    }
  }
  
  cat("\nSummary:\n")
  cat("Total files found:", length(all_files), "\n")
  cat("Successfully processed:", success_count, "\n")
  cat("Failed:", length(all_files) - success_count, "\n")
}

# Command line argument parsing
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
  cat("Usage: Rscript generate_bed_files_from_outliers.R [base_dir] [pattern] [window_size]\n")
  cat("  base_dir: Directory to search for files (default: current directory)\n")
  cat("  pattern: File pattern to match (default: 'overlap_.*\\\\.txt$')\n")
  cat("  window_size: Size of genomic windows in bp (default: 5000)\n")
  quit(status = 1)
}
