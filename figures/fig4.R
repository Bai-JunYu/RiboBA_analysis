#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/main/plot_fig4_core.R",
  run_fun = "run_fig4",
  required = c("dfp", "df_mean"),
  default_rds = "figure_ready_data/fig4_for_plot.rds",
  output_pdf = "figures/Fig_4_MS.pdf"
)
