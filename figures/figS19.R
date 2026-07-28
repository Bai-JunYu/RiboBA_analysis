#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS19_core.R",
  run_fun = "run_fig_s19",
  required = c("panel_A_readlen", "panel_B_hindrance", "panel_C_cutbias",
               "panel_D_add5", "panel_E_ligation", "panel_F_biotype", "panel_G_metagene"),
  default_rds = "figure_ready_data/figS19_for_plot.rds",
  output_pdf = "figures/Fig_S19_riboba_output_examples.pdf"
)
