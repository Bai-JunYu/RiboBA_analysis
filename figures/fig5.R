#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/main/plot_fig5_core.R",
  run_fun = "run_fig5",
  required = c("p2_3", "p_non_3", "p_all"),
  default_rds = "figure_ready_data/fig5_for_plot.rds",
  output_pdf = "figures/Fig_5_fly.pdf"
)
