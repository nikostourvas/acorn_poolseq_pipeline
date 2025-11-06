#!/usr/bin/env python3
"""
Script to extract index column values for each envfactor from CSV file
and save them as separate .txt files.
"""

import csv
import os
from collections import defaultdict

def extract_index_by_envfactor():
    # Input CSV file path
    csv_file = "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/WZA_res/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.01_rho_top2_5.csv"
    
    # Output directory
    output_dir = "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA"
    
    # Dictionary to store index values for each envfactor
    envfactor_indices = defaultdict(list)
    
    # Read the CSV file
    print("Reading CSV file...")
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            envfactor = row['envfactor']
            index_value = row['index']
            envfactor_indices[envfactor].append(index_value)
    
    print(f"Found {len(envfactor_indices)} unique envfactors")
    
    # Process each envfactor
    for envfactor, index_values in envfactor_indices.items():
        print(f"Processing {envfactor}...")
        
        # Create output filename
        output_file = os.path.join(output_dir, f"{envfactor}.txt")
        
        # Write index values to file (one per line)
        with open(output_file, 'w') as f:
            for index_val in index_values:
                f.write(f"{index_val}\n")
        
        print(f"  Saved {len(index_values)} index values to {output_file}")
    
    print("Processing complete!")

if __name__ == "__main__":
    extract_index_by_envfactor()
