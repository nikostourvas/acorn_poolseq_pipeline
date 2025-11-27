library(GenomicRanges)
library(rtracklayer) # for import()
library(data.table)
library(tidyr)
library(dplyr)
library(plyranges)
library(ggplot2)

# Set options to avoid scientific notation in output
# https://divingintogeneticsandgenomics.com/post/three-gotchas-when-using-r-for-genomic-data-analysis/
options(scipen = 500)

# Function to fix formatting issues with LFMM files
fix_lfmm <- function(x) {
  x$SNPid = gsub("Qrob_", "", x$SNPid)
  x$SNPid = gsub("Chr0", "chr", x$SNPid)
  x$SNPid = gsub("Chr", "chr", x$SNPid)
  x$SNPid = gsub("H2.3_", "H2.3.", x$SNPid)
  x = separate(x, col = SNPid, into = c("chr", "pos"), sep = "_")
  x$pos = as.integer(x$pos)
  return(x)
}

# Function to fix formatting issues with BayPass files
fix_baypass <- function(x, envfactor) {
  colnames(x)[1] = "chr"
  colnames(x)[2] = "pos"
  # top 10% of BF values
  # cutoff = x[, quantile(x[[envfactor]], probs = 0.9)]
  # x = x[x[[envfactor]] >= cutoff,]
  # BF threshold
  x = x[x[[envfactor]] >= 10, ]
  x$chr = paste0("chr", x$chr)
  x$chr = gsub("chrQrob_", "", x$chr)
  x$pos = as.integer(x$pos)
  return(x)
}

# Bio1
snps_bio1_lfmm = fread(
  "/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio1/CandidatesOrdered_envbio1_K1_q0.001.csv"
)
snps_bio1_lfmm = fix_lfmm(snps_bio1_lfmm)

snps_bio1_baypass = fread(
  "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate_filtered.csv"
)[, c(1:5, 17)]
snps_bio1_baypass = fix_baypass(snps_bio1_baypass, "bio1_BF")

# Find snps detected by both LFMM and BayPass
snps_bio1_overlap = inner_join(
  snps_bio1_lfmm,
  snps_bio1_baypass,
  by = c("chr", "pos")
)
snp_bio1_gr <- GRanges(
  seqnames = snps_bio1_overlap$chr,
  ranges = IRanges(snps_bio1_overlap$pos, snps_bio1_overlap$pos)
) # 1‑bp ranges

# Load genomic ranges significantly associated with the environment based on LFMM
bed_bio1_lfmm_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio1.bed"
)
snps_bio1_lfmm_gr <- GRanges(
  seqnames = snps_bio1_lfmm$chr,
  ranges = IRanges(snps_bio1_lfmm$pos, snps_bio1_lfmm$pos)
) # 1－bp ranges
hits_bio1_lfmm <- findOverlaps(snps_bio1_lfmm_gr, bed_bio1_lfmm_gr)
snps_bio1_lfmm_in_regions <- snps_bio1_lfmm_gr[queryHits(hits_bio1_lfmm)]

# Load genomic ranges significantly associated with the environment based on BayPass
bed_bio1_baypass_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio1.bed"
)
snps_bio1_baypass_gr <- GRanges(
  seqnames = snps_bio1_baypass$chr,
  ranges = IRanges(snps_bio1_baypass$pos, snps_bio1_baypass$pos)
) # 1－bp ranges
hits_bio1_baypass <- findOverlaps(snps_bio1_baypass_gr, bed_bio1_baypass_gr)
snps_bio1_baypass_in_regions <- snps_bio1_baypass_gr[queryHits(
  hits_bio1_baypass
)]


# Bio2
# snps_bio2_lfmm = fread("/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio2/CandidatesOrdered_envbio2_K1_q0.01.csv")
# snps_bio2_lfmm = fix_lfmm(snps_bio2_lfmm)

# snps_bio2_baypass = fread("/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate.csv")[,c(1:5, 18)]
# snps_bio2_baypass = fix_baypass(snps_bio2_baypass, "bio2_BF")

# # Find snps detected by both LFMM and BayPass
#   snps_bio2_overlap = inner_join(snps_bio2_lfmm, snps_bio2_baypass, by = c("chr", "pos"))
#   snp_bio2_gr <- GRanges(seqnames = snps_bio2_overlap$chr,
#     ranges   = IRanges(snps_bio2_overlap$pos, snps_bio2_overlap$pos))   # 1‑bp ranges

#   # Load genomic ranges significantly associated with the environment based on LFMM
#   bed_bio2_lfmm_gr <- import("/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio2.bed")
#   snps_bio2_lfmm_gr <- GRanges(seqnames = snps_bio2_lfmm$chr,
#     ranges   = IRanges(snps_bio2_lfmm$pos, snps_bio2_lfmm$pos))   # 1－bp ranges
#   hits_bio2_lfmm <- findOverlaps(snps_bio2_lfmm_gr, bed_bio2_lfmm_gr)
#   snps_bio2_lfmm_in_regions <- snps_bio2_lfmm_gr[queryHits(hits_bio2_lfmm)]

#   # Load genomic ranges significantly associated with the environment based on BayPass
#   bed_bio2_baypass_gr <- import("/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio2.bed")
#   snps_bio2_baypass_gr <- GRanges(seqnames = snps_bio2_baypass$chr,
#     ranges   = IRanges(snps_bio2_baypass$pos, snps_bio2_baypass$pos))   # 1－bp ranges
#   hits_bio2_baypass <- findOverlaps(snps_bio2_baypass_gr, bed_bio2_baypass_gr)
#   snps_bio2_baypass_in_regions <- snps_bio2_baypass_gr[queryHits(hits_bio2_baypass)]

# Bio3
snps_bio3_lfmm = fread(
  "/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio3/CandidatesOrdered_envbio3_K1_q0.001.csv"
)
snps_bio3_lfmm = fix_lfmm(snps_bio3_lfmm)
snps_bio3_baypass = fread(
  "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate_filtered.csv"
)[, c(1:5, 19)]
snps_bio3_baypass = fix_baypass(snps_bio3_baypass, "bio3_BF")

# Find snps detected by both LFMM and BayPass
snps_bio3_overlap = inner_join(
  snps_bio3_lfmm,
  snps_bio3_baypass,
  by = c("chr", "pos")
)
snp_bio3_gr <- GRanges(
  seqnames = snps_bio3_overlap$chr,
  ranges = IRanges(snps_bio3_overlap$pos, snps_bio3_overlap$pos)
) # 1‑bp ranges

# Load genomic ranges significantly associated with the environment based on LFMM
bed_bio3_lfmm_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio3.bed"
)
snps_bio3_lfmm_gr <- GRanges(
  seqnames = snps_bio3_lfmm$chr,
  ranges = IRanges(snps_bio3_lfmm$pos, snps_bio3_lfmm$pos)
) # 1－bp ranges
hits_bio3_lfmm <- findOverlaps(snps_bio3_lfmm_gr, bed_bio3_lfmm_gr)
snps_bio3_lfmm_in_regions <- snps_bio3_lfmm_gr[queryHits(hits_bio3_lfmm)]

# Load genomic ranges significantly associated with the environment based on BayPass
bed_bio3_baypass_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio3.bed"
)
snps_bio3_baypass_gr <- GRanges(
  seqnames = snps_bio3_baypass$chr,
  ranges = IRanges(snps_bio3_baypass$pos, snps_bio3_baypass$pos)
) # 1－bp ranges
hits_bio3_baypass <- findOverlaps(snps_bio3_baypass_gr, bed_bio3_baypass_gr)
snps_bio3_baypass_in_regions <- snps_bio3_baypass_gr[queryHits(
  hits_bio3_baypass
)]

# Bio8
# snps_bio8_lfmm = fread("/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio8/CandidatesOrdered_envbio8_K1_q0.01.csv")
# snps_bio8_lfmm = fix_lfmm(snps_bio8_lfmm)

# snps_bio8_baypass = fread("/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate.csv")[,c(1:5, 24)]
# snps_bio8_baypass = fix_baypass(snps_bio8_baypass, "bio8_BF")

# # Find snps detected by both LFMM and BayPass
#   snps_bio8_overlap = inner_join(snps_bio8_lfmm, snps_bio8_baypass, by = c("chr", "pos"))
#   snp_bio8_gr <- GRanges(seqnames = snps_bio8_overlap$chr,
#     ranges   = IRanges(snps_bio8_overlap$pos, snps_bio8_overlap$pos))   # 1‑bp ranges

#   # Load genomic ranges significantly associated with the environment based on LFMM
#   bed_bio8_lfmm_gr <- import("/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio8.bed")
#   snps_bio8_lfmm_gr <- GRanges(seqnames = snps_bio8_lfmm$chr,
#     ranges   = IRanges(snps_bio8_lfmm$pos, snps_bio8_lfmm$pos))   # 1－bp ranges
#   hits_bio8_lfmm <- findOverlaps(snps_bio8_lfmm_gr, bed_bio8_lfmm_gr)
#   snps_bio8_lfmm_in_regions <- snps_bio8_lfmm_gr[queryHits(hits_bio8_lfmm)]

#   # Load genomic ranges significantly associated with the environment based on BayPass
#   bed_bio8_baypass_gr <- import("/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio8.bed")
#   snps_bio8_baypass_gr <- GRanges(seqnames = snps_bio8_baypass$chr,
#     ranges   = IRanges(snps_bio8_baypass$pos, snps_bio8_baypass$pos))   # 1－bp ranges
#   hits_bio8_baypass <- findOverlaps(snps_bio8_baypass_gr, bed_bio8_baypass_gr)
#   snps_bio8_baypass_in_regions <- snps_bio8_baypass_gr[queryHits(hits_bio8_baypass)]

# Bio10
snps_bio10_lfmm = fread(
  "/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio10/CandidatesOrdered_envbio10_K1_q0.01.csv"
)
snps_bio10_lfmm = fix_lfmm(snps_bio10_lfmm)

snps_bio10_baypass = fread(
  "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate_filtered.csv"
)[, c(1:5, 26)]
snps_bio10_baypass = fix_baypass(snps_bio10_baypass, "bio10_BF")

# Find snps detected by both LFMM and BayPass
snps_bio10_overlap = inner_join(
  snps_bio10_lfmm,
  snps_bio10_baypass,
  by = c("chr", "pos")
)
snp_bio10_gr <- GRanges(
  seqnames = snps_bio10_overlap$chr,
  ranges = IRanges(snps_bio10_overlap$pos, snps_bio10_overlap$pos)
) # 1‑bp ranges

# Load genomic ranges significantly associated with the environment based on LFMM
bed_bio10_lfmm_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio10.bed"
)
snps_bio10_lfmm_gr <- GRanges(
  seqnames = snps_bio10_lfmm$chr,
  ranges = IRanges(snps_bio10_lfmm$pos, snps_bio10_lfmm$pos)
) # 1－bp ranges
hits_bio10_lfmm <- findOverlaps(snps_bio10_lfmm_gr, bed_bio10_lfmm_gr)
snps_bio10_lfmm_in_regions <- snps_bio10_lfmm_gr[queryHits(hits_bio10_lfmm)]

# Load genomic ranges significantly associated with the environment based on BayPass
bed_bio10_baypass_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio10.bed"
)
snps_bio10_baypass_gr <- GRanges(
  seqnames = snps_bio10_baypass$chr,
  ranges = IRanges(snps_bio10_baypass$pos, snps_bio10_baypass$pos)
) # 1－bp ranges
hits_bio10_baypass <- findOverlaps(snps_bio10_baypass_gr, bed_bio10_baypass_gr)
snps_bio10_baypass_in_regions <- snps_bio10_baypass_gr[queryHits(
  hits_bio10_baypass
)]

# Bio12
snps_bio12_lfmm = fread(
  "/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio12/CandidatesOrdered_envbio12_K1_q0.001.csv"
)
snps_bio12_lfmm = fix_lfmm(snps_bio12_lfmm)
snps_bio12_baypass = fread(
  "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate_filtered.csv"
)[, c(1:5, 28)]
snps_bio12_baypass = fix_baypass(snps_bio12_baypass, "bio12_BF")
# Find snps detected by both LFMM and BayPass
snps_bio12_overlap = inner_join(
  snps_bio12_lfmm,
  snps_bio12_baypass,
  by = c("chr", "pos")
)
snp_bio12_gr <- GRanges(
  seqnames = snps_bio12_overlap$chr,
  ranges = IRanges(snps_bio12_overlap$pos, snps_bio12_overlap$pos)
) # 1‑bp ranges
# Load genomic ranges significantly associated with the environment based on LFMM
bed_bio12_lfmm_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio12.bed"
)
snps_bio12_lfmm_gr <- GRanges(
  seqnames = snps_bio12_lfmm$chr,
  ranges = IRanges(snps_bio12_lfmm$pos, snps_bio12_lfmm$pos)
) # 1－bp ranges
hits_bio12_lfmm <- findOverlaps(snps_bio12_lfmm_gr, bed_bio12_lfmm_gr)
snps_bio12_lfmm_in_regions <- snps_bio12_lfmm_gr[queryHits(hits_bio12_lfmm)]
# Load genomic ranges significantly associated with the environment based on BayPass
bed_bio12_baypass_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio12.bed"
)
snps_bio12_baypass_gr <- GRanges(
  seqnames = snps_bio12_baypass$chr,
  ranges = IRanges(snps_bio12_baypass$pos, snps_bio12_baypass$pos)
) # 1－bp ranges
hits_bio12_baypass <- findOverlaps(snps_bio12_baypass_gr, bed_bio12_baypass_gr)
snps_bio12_baypass_in_regions <- snps_bio12_baypass_gr[queryHits(
  hits_bio12_baypass
)]


# Bio15
snps_bio15_lfmm = fread(
  "/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio15/CandidatesOrdered_envbio15_K1_q0.001.csv"
)
snps_bio15_lfmm = fix_lfmm(snps_bio15_lfmm)

snps_bio15_baypass = fread(
  "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate_filtered.csv"
)[, c(1:5, 31)]
snps_bio15_baypass = fix_baypass(snps_bio15_baypass, "bio15_BF")

# Find snps detected by both LFMM and BayPass
snps_bio15_overlap = inner_join(
  snps_bio15_lfmm,
  snps_bio15_baypass,
  by = c("chr", "pos")
)
snp_bio15_gr <- GRanges(
  seqnames = snps_bio15_overlap$chr,
  ranges = IRanges(snps_bio15_overlap$pos, snps_bio15_overlap$pos)
) # 1‑bp ranges

# Load genomic ranges significantly associated with the environment based on LFMM
bed_bio15_lfmm_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio15.bed"
)
snps_bio15_lfmm_gr <- GRanges(
  seqnames = snps_bio15_lfmm$chr,
  ranges = IRanges(snps_bio15_lfmm$pos, snps_bio15_lfmm$pos)
) # 1－bp ranges
hits_bio15_lfmm <- findOverlaps(snps_bio15_lfmm_gr, bed_bio15_lfmm_gr)
snps_bio15_lfmm_in_regions <- snps_bio15_lfmm_gr[queryHits(hits_bio15_lfmm)]

# Load genomic ranges significantly associated with the environment based on BayPass
bed_bio15_baypass_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio15.bed"
)
snps_bio15_baypass_gr <- GRanges(
  seqnames = snps_bio15_baypass$chr,
  ranges = IRanges(snps_bio15_baypass$pos, snps_bio15_baypass$pos)
) # 1－bp ranges
hits_bio15_baypass <- findOverlaps(snps_bio15_baypass_gr, bed_bio15_baypass_gr)
snps_bio15_baypass_in_regions <- snps_bio15_baypass_gr[queryHits(
  hits_bio15_baypass
)]

# Bio16
# snps_bio16_lfmm = fread(
#   "/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio16/CandidatesOrdered_envbio16_K1_q0.01.csv"
# )
# snps_bio16_lfmm = fix_lfmm(snps_bio16_lfmm)

# snps_bio16_baypass = fread(
#   "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate.csv"
# )[, c(1:5, 32)]
# snps_bio16_baypass = fix_baypass(snps_bio16_baypass, "bio16_BF")

# # Find snps detected by both LFMM and BayPass
# snps_bio16_overlap = inner_join(
#   snps_bio16_lfmm,
#   snps_bio16_baypass,
#   by = c("chr", "pos")
# )
# snp_bio16_gr <- GRanges(
#   seqnames = snps_bio16_overlap$chr,
#   ranges = IRanges(snps_bio16_overlap$pos, snps_bio16_overlap$pos)
# ) # 1‑bp ranges

# # Load genomic ranges significantly associated with the environment based on LFMM
# bed_bio16_lfmm_gr <- import(
#   "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio16.bed"
# )
# snps_bio16_lfmm_gr <- GRanges(
#   seqnames = snps_bio16_lfmm$chr,
#   ranges = IRanges(snps_bio16_lfmm$pos, snps_bio16_lfmm$pos)
# ) # 1－bp ranges
# hits_bio16_lfmm <- findOverlaps(snps_bio16_lfmm_gr, bed_bio16_lfmm_gr)
# snps_bio16_lfmm_in_regions <- snps_bio16_lfmm_gr[queryHits(hits_bio16_lfmm)]

# # Load genomic ranges significantly associated with the environment based on BayPass
# bed_bio16_baypass_gr <- import(
#   "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio16.bed"
# )
# snps_bio16_baypass_gr <- GRanges(
#   seqnames = snps_bio16_baypass$chr,
#   ranges = IRanges(snps_bio16_baypass$pos, snps_bio16_baypass$pos)
# ) # 1－bp ranges
# hits_bio16_baypass <- findOverlaps(snps_bio16_baypass_gr, bed_bio16_baypass_gr)
# snps_bio16_baypass_in_regions <- snps_bio16_baypass_gr[queryHits(
#   hits_bio16_baypass
# )]

# Bio18
snps_bio18_lfmm = fread(
  "/home/geneticsShare/LFMM_ANALYSES/res_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot/Full_analysis_K_1/environment_bio18/CandidatesOrdered_envbio18_K1_q0.01.csv"
)
snps_bio18_lfmm = fix_lfmm(snps_bio18_lfmm)
snps_bio18_baypass = fread(
  "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_concatenated_res_covariate_filtered.csv"
)[, c(1:5, 34)]
snps_bio18_baypass = fix_baypass(snps_bio18_baypass, "bio18_BF")
# Find snps detected by both LFMM and BayPass
snps_bio18_overlap = inner_join(
  snps_bio18_lfmm,
  snps_bio18_baypass,
  by = c("chr", "pos")
)
snp_bio18_gr <- GRanges(
  seqnames = snps_bio18_overlap$chr,
  ranges = IRanges(snps_bio18_overlap$pos, snps_bio18_overlap$pos)
) # 1‑bp ranges
# Load genomic ranges significantly associated with the environment based on LFMM
bed_bio18_lfmm_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_LFMM/q0.01_1CA/bio18.bed"
)
snps_bio18_lfmm_gr <- GRanges(
  seqnames = snps_bio18_lfmm$chr,
  ranges = IRanges(snps_bio18_lfmm$pos, snps_bio18_lfmm$pos)
) # 1－bp ranges
hits_bio18_lfmm <- findOverlaps(snps_bio18_lfmm_gr, bed_bio18_lfmm_gr)
snps_bio18_lfmm_in_regions <- snps_bio18_lfmm_gr[queryHits(hits_bio18_lfmm)]
# Load genomic ranges significantly associated with the environment based on BayPass
bed_bio18_baypass_gr <- import(
  "/home/geneticsShare/results/comparisons_wp3_wp4/Robust_BayPass/q0.01_rho_top_2_5_1CA/bio18.bed"
)
snps_bio18_baypass_gr <- GRanges(
  seqnames = snps_bio18_baypass$chr,
  ranges = IRanges(snps_bio18_baypass$pos, snps_bio18_baypass$pos)
) # 1－bp ranges
hits_bio18_baypass <- findOverlaps(snps_bio18_baypass_gr, bed_bio18_baypass_gr)
snps_bio18_baypass_in_regions <- snps_bio18_baypass_gr[queryHits(
  hits_bio18_baypass
)]


# Bring in mean MAF values to further filter the SNPs
maf_gr = fread(
  "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/data/WZA/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_WZA_input.csv"
)[, c(1:2, 67)] |>
  rename(chr = CHR, pos = POS) |>
  mutate(chr = paste0("chr", chr)) |>
  mutate(chr = gsub("chrQrob_", "", chr)) |>
  mutate(seqnames = chr, start = pos, end = pos) |>
  as_granges()

# Add metadata to each GRanges object before merging
snps_bio1_lfmm_in_regions$bioclim <- "bio1"
snps_bio1_lfmm_in_regions$method <- "lfmm"

snps_bio1_baypass_in_regions$bioclim <- "bio1"
snps_bio1_baypass_in_regions$method <- "baypass"

snps_bio3_lfmm_in_regions$bioclim <- "bio3"
snps_bio3_lfmm_in_regions$method <- "lfmm"

snps_bio3_baypass_in_regions$bioclim <- "bio3"
snps_bio3_baypass_in_regions$method <- "baypass"

snps_bio10_lfmm_in_regions$bioclim <- "bio10"
snps_bio10_lfmm_in_regions$method <- "lfmm"

snps_bio10_baypass_in_regions$bioclim <- "bio10"
snps_bio10_baypass_in_regions$method <- "baypass"

snps_bio12_lfmm_in_regions$bioclim <- "bio12"
snps_bio12_lfmm_in_regions$method <- "lfmm"

snps_bio12_baypass_in_regions$bioclim <- "bio12"
snps_bio12_baypass_in_regions$method <- "baypass"

snps_bio15_lfmm_in_regions$bioclim <- "bio15"
snps_bio15_lfmm_in_regions$method <- "lfmm"

snps_bio15_baypass_in_regions$bioclim <- "bio15"
snps_bio15_baypass_in_regions$method <- "baypass"

snps_bio18_lfmm_in_regions$bioclim <- "bio18"
snps_bio18_lfmm_in_regions$method <- "lfmm"

snps_bio18_baypass_in_regions$bioclim <- "bio18"
snps_bio18_baypass_in_regions$method <- "baypass"

# Merge genomic ranges for different envfactors
# merged_gr = c(
#   snps_bio1_lfmm_in_regions,
#   snps_bio1_baypass_in_regions,
#   snps_bio3_lfmm_in_regions,
#   snps_bio3_baypass_in_regions,
#   snps_bio10_lfmm_in_regions,
#   snps_bio10_baypass_in_regions,
#   snps_bio12_lfmm_in_regions,
#   snps_bio12_baypass_in_regions,
#   snps_bio15_lfmm_in_regions,
#   snps_bio15_baypass_in_regions,
#   snps_bio18_lfmm_in_regions,
#   snps_bio18_baypass_in_regions
# ) |>
#   unique()

# After adding metadata, combine methods for duplicate SNPs
merge_ranges <- function(gr) {
  gr |>
    as_tibble() |>
    group_by(seqnames, start) |>
    summarize(
      end = first(end),
      bioclim = paste(unique(bioclim), collapse = ","),
      method = paste(unique(method), collapse = ","),
      .groups = "drop"
    ) |>
    as_granges()
}
merged_gr <- merge_ranges(c(
  snps_bio1_lfmm_in_regions,
  snps_bio1_baypass_in_regions,
  snps_bio3_lfmm_in_regions,
  snps_bio3_baypass_in_regions,
  snps_bio10_lfmm_in_regions,
  snps_bio10_baypass_in_regions,
  snps_bio12_lfmm_in_regions,
  snps_bio12_baypass_in_regions,
  snps_bio15_lfmm_in_regions,
  snps_bio15_baypass_in_regions,
  snps_bio18_lfmm_in_regions,
  snps_bio18_baypass_in_regions
))

# Join with MAF values and filter for mean MAF >= 0.05
merged_gr <- merged_gr |>
  join_overlap_intersect(maf_gr)

filtered_gr <- merged_gr |>
  filter(MAF >= 0.45)

filtered_gr

# Further filter: Keep the highest-MAF SNP in each 1000-bp window
# gr: a GRanges of single‑base SNPs
# min_dist: minimum distance (in bp) between kept SNPs
# Assumes mcols(gr)$MAF exists
thin_snps_by_maxMAF <- function(gr, min_dist) {
  # 1) Build blocks that merge any two SNPs with gap < min_dist
  blocks <- GenomicRanges::reduce(gr, min.gapwidth = min_dist)

  # 2) Find which SNPs fall into which block
  hits <- findOverlaps(blocks, gr)
  snps_in_block <- split(subjectHits(hits), queryHits(hits))

  # 3) In each block, pick the SNP with highest MAF
  maf <- mcols(gr)$MAF
  best_per_block <- sapply(snps_in_block, function(idxs) {
    # which.max returns the first max if ties
    idxs[which.max(maf[idxs])]
  })

  # 4) Return those SNPs, sorted by genomic position
  sort(gr[best_per_block])
}

# keep SNPs ≥x bp apart, picking highest‐MAF in each cluster
thin_hiMAF_snps <- thin_snps_by_maxMAF(filtered_gr, min_dist = 5000)
print(as_tibble(thin_hiMAF_snps), n = 1000)


######### Apply thin function for each bioclimatic variable separately ########
# merge snps_bio1_lfmm_gr, snps_bio1_baypass_gr
thin_hiMAF_snps_bio1 <- merge_ranges(
  c(snps_bio1_lfmm_in_regions, snps_bio1_baypass_in_regions)
) |>
  join_overlap_intersect(maf_gr) |>
  filter(MAF >= 0.25) |>
  thin_snps_by_maxMAF(min_dist = 5000)
print(as_tibble(thin_hiMAF_snps_bio1), n = 1000)

# merge snps_bio3_lfmm_gr, snps_bio3_baypass_gr
thin_hiMAF_snps_bio3 <- merge_ranges(
  c(snps_bio3_lfmm_in_regions, snps_bio3_baypass_in_regions)
) |>
  join_overlap_intersect(maf_gr) |>
  thin_snps_by_maxMAF(min_dist = 5000) |>
  filter(MAF >= 0.25)
print(as_tibble(thin_hiMAF_snps_bio3), n = 1000)

# merge snps_bio10_lfmm_gr, snps_bio10_baypass_gr
thin_hiMAF_snps_bio10 <- merge_ranges(
  c(snps_bio10_lfmm_in_regions, snps_bio10_baypass_in_regions)
) |>
  join_overlap_intersect(maf_gr) |>
  thin_snps_by_maxMAF(min_dist = 5000) |>
  filter(MAF >= 0.25)
print(as_tibble(thin_hiMAF_snps_bio10), n = 1000)

# merge snps_bio12_lfmm_gr, snps_bio12_baypass_gr
thin_hiMAF_snps_bio12 <- merge_ranges(
  c(snps_bio12_lfmm_in_regions, snps_bio12_baypass_in_regions)
) |>
  join_overlap_intersect(maf_gr) |>
  thin_snps_by_maxMAF(min_dist = 5000) |>
  filter(MAF >= 0.25)
print(as_tibble(thin_hiMAF_snps_bio12), n = 1000)

# merge snps_bio15_lfmm_gr, snps_bio15_baypass_gr
thin_hiMAF_snps_bio15 <- merge_ranges(
  c(snps_bio15_lfmm_in_regions, snps_bio15_baypass_in_regions)
) |>
  join_overlap_intersect(maf_gr) |>
  thin_snps_by_maxMAF(min_dist = 5000) |>
  filter(MAF >= 0.25)
print(as_tibble(thin_hiMAF_snps_bio15), n = 1000)

# merge snps_bio18_lfmm_gr, snps_bio18_baypass_gr
thin_hiMAF_snps_bio18 <- merge_ranges(
  c(snps_bio18_lfmm_in_regions, snps_bio18_baypass_in_regions)
) |>
  join_overlap_intersect(maf_gr) |>
  thin_snps_by_maxMAF(min_dist = 5000) |>
  filter(MAF >= 0.25)
print(as_tibble(thin_hiMAF_snps_bio18), n = 1000)


# Combine the individually filtered datasets
combined_ranges <- merge_ranges(
  c(
    thin_hiMAF_snps_bio1,
    thin_hiMAF_snps_bio3,
    thin_hiMAF_snps_bio10,
    thin_hiMAF_snps_bio12,
    thin_hiMAF_snps_bio15,
    thin_hiMAF_snps_bio18
  )
)

print(as_tibble(combined_ranges), n = 1000)

# Thin the combined dataset. No SNPs should be closer than 5,000 bp.
# MAF is NOT used here.
combined_test = combined_ranges |>
  join_overlap_intersect(maf_gr) |>
  thin_snps_by_maxMAF(min_dist = 500000)
print(as_tibble(combined_test), n = 1000)

# Bio12 and Bio15 and possibly Bio3 still have too many SNPs. Filter only these rows using MAF. The rest of the rows should remain intact.
combined_balanced = combined_test |>
  as_tibble() |>
  mutate(bioclim = as.character(bioclim)) |>
  group_by(bioclim) |>
  mutate(
    MAF_threshold = case_when(
      bioclim == "bio12" ~ 0.25,
      bioclim == "bio15" ~ 0.25,
      TRUE ~ 0.25
    )
  ) |>
  filter(MAF >= MAF_threshold) |>
  select(-MAF_threshold) |>
  ungroup() |>
  as_granges()


print(as_tibble(combined_balanced), n = 1000)


# Still too many SNPs. Subsample them so that there are 26 SNPs per bioclimatic variable
set.seed(42)
thin_balanced_snps <- combined_balanced |>
  as_tibble() |>
  group_by(bioclim) |>
  slice_sample(prop = 1) |> # shuffle within each group
  slice_head(n = 22) |> # take up to 30 from each group
  ungroup() |>
  as_granges()
print(as_tibble(thin_balanced_snps), n = 1000)
##################################################333

# Export results from original method
as_tibble(thin_hiMAF_snps) |>
  transmute(
    chr = as.character(seqnames),
    pos = start,
    # rsid = if_else(is.na(ID), ".", ID),
    # ref  = REF,
    # alt  = ALT,
    MAF,
    bioclim,
    method
  ) |>
  fwrite(
    "2025-08-13_massarray_proposed_snps_bio1-3-10-12-15-18_q0001-bf10_minMAF-038_original.tsv",
    sep = "\t"
  )

# Export results from balanced method
as_tibble(thin_balanced_snps) |>
  transmute(
    chr = as.character(seqnames),
    pos = start,
    # rsid = if_else(is.na(ID), ".", ID),
    # ref  = REF,
    # alt  = ALT,
    MAF,
    bioclim,
    method
  ) |>
  fwrite(
    "2025-08-13_massarray_proposed_snps_bio1-3-10-12-15-18_q0001-bf10_minMAF-038_balanced.tsv",
    sep = "\t"
  )

# export unfiltered merged results for reference
as_tibble(merged_gr) |>
  transmute(
    chr = as.character(seqnames),
    pos = start,
    # rsid = if_else(is.na(ID), ".", ID),
    # ref  = REF,
    # alt  = ALT,
    MAF,
    bioclim,
    method
  ) |>
  fwrite(
    "2025-08-13_massarray_proposed_snps_bio1-3-10-12-15-18_q0001-bf10_unfiltered.tsv",
    sep = "\t"
  )

# export combined_test results for reference
as_tibble(combined_balanced) |>
  transmute(
    chr = as.character(seqnames),
    pos = start,
    # rsid = if_else(is.na(ID), ".", ID),
    # ref  = REF,
    # alt  = ALT,
    MAF,
    bioclim,
    method
  ) |>
  fwrite(
    "2025-08-25_massarray_proposed_snps_bio1-3-10-12-15-18_q0001-bf10_variable-minMAF.tsv",
    sep = "\t"
  )

as_tibble(thin_balanced_snps) |>
  transmute(
    chr = as.character(seqnames),
    pos = start,
    # rsid = if_else(is.na(ID), ".", ID),
    # ref  = REF,
    # alt  = ALT,
    MAF,
    bioclim,
    method
  ) |>
  fwrite(
    "2025-11-27_massarray_proposed_snps_bio1-3-10-12-15-18_q001-bf10_balanced-26per-variable.tsv",
    sep = "\t"
  )

# Explore how many SNPs I have per bioclimatic variable
# Count the number of SNPs per bioclimatic variable
count_per_bioclim <- as_tibble(thin_hiMAF_snps) |>
  group_by(bioclim) |>
  summarise(n = n(), .groups = "drop")

# make a bar plot
ggplot(count_per_bioclim, aes(x = bioclim, y = n)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Number of SNPs per Bioclimatic Variable",
    x = "Bioclimatic Variable",
    y = "Number of SNPs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# same plot for the merged_gr dataset
count_per_bioclim_merged <- as_tibble(merged_gr) |>
  group_by(bioclim) |>
  summarise(n = n(), .groups = "drop")

# add a new column to indicate whether one or multple bioclim variables are present. Write "multi" if multiple bioclim variables are present, otherwise write "single" in the new column
count_per_bioclim_merged$bioclim_count <- ifelse(
  grepl(",", count_per_bioclim_merged$bioclim),
  "Multiple Env Factors\nfor the same SNP",
  "Single Env Factor"
)

ggplot(count_per_bioclim_merged, aes(x = bioclim, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap(~bioclim_count, scales = "free_x") +
  labs(
    title = "Number of SNPs per Bioclimatic Variable (Unfiltered)",
    x = "Bioclimatic Variable",
    y = "Number of SNPs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# same plot for combined_test dataset
count_per_bioclim_combined_test <- as_tibble(combined_balanced) |>
  group_by(bioclim) |>
  summarise(n = n(), .groups = "drop")

# add a new column to indicate whether one or multple bioclim variables are present. Write "multi" if multiple bioclim variables are present, otherwise write "single" in the new column
count_per_bioclim_combined_test$bioclim_count <- ifelse(
  grepl(",", count_per_bioclim_combined_test$bioclim),
  "Multiple Env Factors\nfor the same SNP",
  "Single Env Factor"
)

count_per_bioclim_combined_test$minMAF <- ifelse(
  count_per_bioclim_combined_test$bioclim == "bio12",
  0.33,
  ifelse(
    count_per_bioclim_combined_test$bioclim == "bio15",
    0.45,
    0.20
  )
)

ggplot(
  count_per_bioclim_combined_test,
  aes(x = bioclim, y = n, fill = as.character(minMAF))
) +
  geom_bar(stat = "identity") +
  facet_wrap(~bioclim_count, scales = "free_x") +
  scale_fill_manual(values = c("red", "orange", "lightblue")) +
  labs(
    title = "Number of SNPs per Bioclimatic Variable",
    x = "Bioclimatic Variable",
    y = "Number of SNPs",
    fill = "Minimum\nMAF "
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# same plot for thin_balanced_snps dataset
count_per_bioclim_thin_balanced <- as_tibble(thin_balanced_snps) |>
  group_by(bioclim) |>
  summarise(n = n(), .groups = "drop")

# add a new column to indicate whether one or multple bioclim variables are present. Write "multi" if multiple bioclim variables are present, otherwise write "single" in the new column
count_per_bioclim_thin_balanced$bioclim_count <- ifelse(
  grepl(",", count_per_bioclim_thin_balanced$bioclim),
  "Multiple Env Factors\nfor the same SNP",
  "Single Env Factor"
)

ggplot(
  count_per_bioclim_thin_balanced,
  aes(x = bioclim, y = n)
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  facet_wrap(~bioclim_count, scales = "free_x") +
  labs(
    title = "Number of SNPs per Bioclimatic Variable (Balanced)",
    x = "Bioclimatic Variable",
    y = "Number of SNPs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "2025-11-06_bioclim_variable_snp_counts_barplot.png",
  width = 10,
  height = 6,
  dpi = 300
)

# Make a manhattan plot of the filtered SNPs and another one for the unfiltered SNPs
thin_hiMAF_snps_df <- as_tibble(thin_hiMAF_snps) |>
  mutate(pos = start) |>
  select(seqnames, pos, MAF, bioclim, method)
thin_hiMAF_snps_df$seqnames <- factor(
  thin_hiMAF_snps_df$seqnames,
  levels = paste0("chr", 1:12)
)

merged_gr_df <- as_tibble(merged_gr) |>
  mutate(pos = start) |>
  select(seqnames, pos, MAF, bioclim, method)
merged_gr_df$seqnames <- factor(
  merged_gr_df$seqnames,
  levels = paste0("chr", 1:12)
)

create_manhattan_plot <- function(data, title) {
  ggplot(data, aes(x = pos, y = MAF)) +
    geom_point(aes(color = bioclim), alpha = 0.7) +
    facet_wrap(~seqnames, scales = "free_x", nrow = 1) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = title,
      x = "Position",
      y = "MAF"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}
# Create Manhattan plots for filtered and unfiltered datasets
manhattan_filtered <- create_manhattan_plot(
  thin_hiMAF_snps_df,
  "Manhattan Plot of Filtered SNPs (MAF >= 0.38)"
)
manhattan_filtered

manhattan_unfiltered <- create_manhattan_plot(
  merged_gr_df,
  "Manhattan Plot of Unfiltered SNPs"
)
manhattan_unfiltered
