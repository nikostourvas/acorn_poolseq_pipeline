setwd(
    "/home/geneticsShare/results/comparisons_wp3_wp4/Paired_GEA_overlap_with_Robust_Models/"
)

#1CA
# example_results <- compare_three_tables(
#     table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.1_rho_top2_5.csv",
#     table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv",
#     table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.1.csv",
#     table1_name = "BayPass",
#     table2_name = "LFMM",
#     table3_name = "Paired_GEA",
#     output_dir = "1CA_q0.1",
#     n_windows = 137408,
#     generate_venn_diagram = TRUE)

# example_results <- compare_three_tables(
#     table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.001_rho_top2_5.csv",
#     table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.001.csv",
#     table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.001.csv",
#     table1_name = "BayPass",
#     table2_name = "LFMM",
#     table3_name = "Paired_GEA",
#     output_dir = "1CA_q0.001",
#     n_windows = 137408,
#     generate_venn_diagram = TRUE)

example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/WZA_res/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.1_noGIF_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/res_robust/WZA_res/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1_noGIF.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_rmTEs_significant_windows_pairedGEA_q0.1_noGIF.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "1CA_q0.1_noGIF",
    n_windows = 104874,
    generate_venn_diagram = TRUE
)

example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/BayPassACORN_AUTh/results/WZA_res/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_q0.001_noGIF_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/res_robust/WZA_res/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.001_noGIF.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpet_1CA_MinDP20_MaxMeanDP178_Miss0075_MAF005_HDplot_rmTEs_significant_windows_pairedGEA_q0.001_noGIF.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "1CA_q0.001_noGIF",
    n_windows = 104874,
    generate_venn_diagram = TRUE
)

#1EA
example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_q0.1_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.1.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "1EA_q0.1",
    n_windows = 131680,
    generate_venn_diagram = TRUE
)

example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_q0.001_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.001.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpet_1EA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.001.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "1EA_q0.001",
    n_windows = 131680,
    generate_venn_diagram = TRUE
)

#2CA
example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_q0.1_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpub_2CP_MinDP20_MaxMeanDP166_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.1.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "2CA_q0.1",
    n_windows = 130043,
    generate_venn_diagram = TRUE
)

example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_q0.001_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2CA_MinDP20_MaxMeanDP165_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.001.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpub_2CP_MinDP20_MaxMeanDP166_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.001.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "2CA_q0.001",
    n_windows = 130043,
    generate_venn_diagram = TRUE
)

#2EA
example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_q0.1_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpub_2EP_MinDP20_MaxMeanDP175_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.1.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "2EA_q0.1",
    n_windows = 128696,
    generate_venn_diagram = TRUE
)

example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_q0.001_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qpub_2EA_MinDP20_MaxMeanDP171_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.001.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qpub_2EP_MinDP20_MaxMeanDP175_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.001.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "2EA_q0.001",
    n_windows = 128696,
    generate_venn_diagram = TRUE
)

#3CA
example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_q0.1_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.1.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "3CA_q0.1",
    n_windows = 130426,
    generate_venn_diagram = TRUE
)

example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_q0.001_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.001.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qrob_3CA_MinDP20_MaxMeanDP189_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.001.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "3CA_q0.001",
    n_windows = 130426,
    generate_venn_diagram = TRUE
)

#3EA
example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_q0.1_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.1.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qrob_3EP_MinDP20_MaxMeanDP196_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.1.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "3EA_q0.1",
    n_windows = 130326,
    generate_venn_diagram = TRUE
)

example_results <- compare_three_tables(
    table1_path = "/home/geneticsShare/BayPassAcorn/parametric_robust/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_q0.001_rho_top2_5.csv",
    table2_path = "/home/geneticsShare/LFMM_ANALYSES/parametric_robust/ACORN_Qrob_3EA_MinDP20_MaxMeanDP191_Miss0075_MAF005_HDplot_significant_windows_LFMM_q0.001.csv",
    table3_path = "/home/geneticsShare/ACORN_paired_GEA/results/WZA/ACORN_Qrob_3EP_MinDP20_MaxMeanDP196_Miss0075_MAF005_HDplot_significant_windows_pairedGEA_q0.001.csv",
    table1_name = "BayPass",
    table2_name = "LFMM",
    table3_name = "Paired_GEA",
    output_dir = "3EA_q0.001",
    n_windows = 130326,
    generate_venn_diagram = TRUE
)
