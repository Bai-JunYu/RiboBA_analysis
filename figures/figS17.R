#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS17_core.R",
  run_fun = "run_fig_s17",
  required = c("deltas", "tests"),
  default_rds = "figure_ready_data/figS17_for_plot.rds",
  output_pdf = "figures/Fig_S17_displayed_grid_auc_paired_tests.pdf"
)
