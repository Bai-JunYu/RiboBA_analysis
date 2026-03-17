#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS4_core.R",
  run_fun = "run_fig_s4",
  required = character(0),
  default_rds = NULL,
  output_pdf = "figures/Fig_S4_runtime.pdf"
)
