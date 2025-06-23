# ── Load the required libraries ───────────────────────────────────────────
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(VennDiagram)
library(data.table)

# ── USER CONFIGURATION ────────────────────────────────────────────────────
# Set the desired q-value threshold (this will be used in file paths and output names)
Q_THRESHOLD <- "q0.1"  # Change this value as needed (e.g., "q0.05", "q0.01")

my_pal <- c(
  neutral_light = "#BBBBBB", neutral_dark = "#000000",
  blue_light    = "#6699CC", blue_dark    = "#004488",
  gold_light    = "#EECC66", gold_dark    = "#997700",
  sig_red       = "#CC3311"
)

# functions
## helper to read a single file and add a dataset ID column ---------------
read_one <- function(f) {
  fread(f, select = c("index", "envfactor"))[
    , dataset := tools::file_path_sans_ext(basename(f))]
}

pairwise_overlaps_vec <- function(list_of_vectors, n_windows) {
  combn(names(list_of_vectors), 2, simplify = FALSE) |>
    map_dfr(\(pair) {
      a <- list_of_vectors[[pair[1]]]
      b <- list_of_vectors[[pair[2]]]
      tibble(
        set_A   = pair[1],
        set_B   = pair[2],
        n_A     = length(a),
        n_B     = length(b),
        n_AB    = length(intersect(a, b)),
        prop_AB = n_AB / min(n_A, n_B),
        pvalue  = phyper(q = length(intersect(a, b)),
                         m = length(a),
                         n = n_windows - length(a),
                         k = length(b),
                         lower.tail=FALSE)
      )
    })
}

find_overlap <- function(files, dataset, n_windows){
        all_dt <- rbindlist(map(files, read_one))   # files = files_1 ∪ files_2
        
        # fix naming: replace "total.solar.radiation.method.1" with "total.solar.radiation" in all_dt$envfactor
        all_dt[envfactor == "total.solar.radiation.method.1", envfactor := "total.solar.radiation"]
        # remove any rows with envfactor == "diffus.solar.radiation" or envfactor == "direct.solar.radiation"
        all_dt <- all_dt[!envfactor %in% c("diffus.solar.radiation", "direct.solar.radiation")]
        
        results <- all_dt[
          , .(index_list = list(unique(index))),      # vector per dataset
          by = .(envfactor, dataset)                 # grouped by both
        ][
          , {                                         # second data.table step
            # turn each envfactor group into a named list of vectors
            idx_sets <- set_names(index_list, dataset)
            pairwise_overlaps_vec(idx_sets, n_windows)         # run the overlap function
          },
          by = envfactor
        ]
        
        # write one CSV per envfactor ─────────────────────────────────
        ## create an output folder
        out_dir <- paste0("pairwise_overlap_", Q_THRESHOLD, "_", dataset)
        dir.create(out_dir, showWarnings = FALSE)
        
        walk(unique(results$envfactor), \(ef) {
          filename <- str_c(dataset, "_", str_replace_all(ef, "[^A-Za-z0-9]+", "_"), "_", Q_THRESHOLD, "_pairwise_overlap.csv")
          results |> filter(envfactor == ef) |>
            write_csv(file = file.path(out_dir, filename))
        })
        
        # all_dt  already has:  dataset | envfactor | index
        # ---------------------------------------------------------------------------
        library(data.table)
        library(stringr)
        
        # Build a table that holds only the *intersection* per envfactor ───────
        common_by_env <- all_dt[
          , .(idx = list(unique(index))),           # one vector per (dataset × envfactor)
          by = .(envfactor, dataset)
        ][
          , {                                       # now inside a single envfactor
            common_idx <- if (.N > 1)             # ≥ 2 datasets -> real intersection
              Reduce(intersect, idx)
            else                     # only one dataset -> keep all indices
              unlist(idx)
            .(index = common_idx)                 # data.table will unnest into rows
          },
          by = envfactor
        ]
        
        # Write one file per envfactor ─────────────────────────────────────────
        ## create an output folder
        out_dir <- paste0("common_genomic_windows_", Q_THRESHOLD, "_", dataset)
        dir.create(out_dir, showWarnings = FALSE)
        
        common_by_env[
          , {
            filename <- str_c(dataset, "_", str_replace_all(envfactor, '[^A-Za-z0-9]+', '_'), "_", Q_THRESHOLD, "_common_index.txt")
            fname <- file.path(out_dir, filename)
            fwrite(.SD[, .(index)], fname, col.names = FALSE)
          },
          by = envfactor
        ]
        
        message("Common‑index files written to ", out_dir, "/")
}

# ── Run the function for all datasets ─────────────────────────────────────
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "1CA"
find_overlap(files, dataset, 137408)
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "1EA"
find_overlap(files, dataset, 131680)
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpet_1GA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1GA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "1GA"
find_overlap(files, dataset, 133866)
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "2CA"
find_overlap(files, dataset, 130043)
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "2EA"
find_overlap(files, dataset, 128696)
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpub_2GA_MinDP20_MaxMeanDP164_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2GA_MinDP20_MaxMeanDP164_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "2GA"
find_overlap(files, dataset, 122034)
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "3CA"
find_overlap(files, dataset, 130426)
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "3EA"
find_overlap(files, dataset, 130321)
files <- c(paste0("/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qrob_3GA_MinDP20_MaxMeanDP188_Miss0075_MAF005_HDplot_significant_windows_", Q_THRESHOLD, "_rho_top2_5.csv"),
           paste0("/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3GA_MinDP20_MaxMeanDP188_Miss0075_MAF005_HDplot_significant_windows_LFMM_", Q_THRESHOLD, ".csv"))
dataset <- "3GA"
find_overlap(files, dataset, 129258)





# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "1CA"
# find_overlap(files, dataset, 137408)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "1EA"
# find_overlap(files, dataset, 131680)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpet_1GA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qpet_1GA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "1GA"
# find_overlap(files, dataset, 133866)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "2CA"
# find_overlap(files, dataset, 130043)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "2EA"
# find_overlap(files, dataset, 128696)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpub_2GA_MinDP20_MaxMeanDP164_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qpub_2GA_MinDP20_MaxMeanDP164_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "2GA"
# find_overlap(files, dataset, 122034)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "3CA"
# find_overlap(files, dataset, 130426)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "3EA"
# find_overlap(files, dataset, 130321)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qrob_3GA_MinDP20_MaxMeanDP188_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric/ACORN_Qrob_3GA_MinDP20_MaxMeanDP188_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "3GA"
# find_overlap(files, dataset, 129258)










# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "1CA"
# find_overlap(files, dataset, 137408)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "1EA"
# find_overlap(files, dataset, 131680)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpet_1GA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1GA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "1GA"
# find_overlap(files, dataset, 133866)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "2CA"
# find_overlap(files, dataset, 130043)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "2EA"
# find_overlap(files, dataset, 128696)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qpub_2GA_MinDP20_MaxMeanDP164_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2GA_MinDP20_MaxMeanDP164_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "2GA"
# find_overlap(files, dataset, 122034)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "3CA"
# find_overlap(files, dataset, 130426)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "3EA"
# find_overlap(files, dataset, 130321)
# files <- c("/home/geneticsShare/BayPassAcorn/parametric/ACORN_Qrob_3GA_MinDP20_MaxMeanDP188_Miss0075_MAF005_HDplot_significant_windows_q0.1.csv",
#            "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3GA_MinDP20_MaxMeanDP188_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv")
# dataset <- "3GA"
# find_overlap(files, dataset, 129258)
