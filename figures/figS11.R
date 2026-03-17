#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS11_core.R",
  run_fun = "run_fig_s11",
  required = c("res1", "res2", "res3", "res1o", "res2o", "res3o"),
  default_rds = "figure_ready_data/figS11_for_plot.rds",
  output_pdf = "figures/Fig_S11_fly_scatter_dot.pdf"
)
