#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS12_core.R",
  run_fun = "run_fig_s12",
  required = c("meta_profile", "offset_exact_accuracy"),
  default_rds = "figure_ready_data/figS12_for_plot.rds",
  output_pdf = "figures/Fig_S12_simulated_psite_assignment_periodicity.pdf"
)
