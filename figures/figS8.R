#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS8_core.R",
  run_fun = "run_fig_s8",
  required = c("df_ms_cnt"),
  default_rds = "figure_ready_data/figS8_for_plot.rds",
  output_pdf = "figures/Fig_S8_MS_number.pdf"
)
