#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS6_core.R",
  run_fun = "run_fig_s6",
  required = c("p_bias", "p_corr", "plig"),
  default_rds = "figure_ready_data/figS6_for_plot.rds",
  output_pdf = "figures/Fig_S6_bias_corr.pdf"
)
