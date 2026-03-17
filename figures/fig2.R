#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/main/plot_fig2_core.R",
  run_fun = "run_fig2",
  required = c("p_roc1", "p_pr", "p_bar1"),
  default_rds = "figure_ready_data/fig2_for_plot.rds",
  output_pdf = "figures/Fig_2_norf_simulation.pdf"
)
