#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS7_core.R",
  run_fun = "run_fig_s7",
  required = c("df_cnt"),
  default_rds = "figure_ready_data/figS7_for_plot.rds",
  output_pdf = "figures/Fig_S7_9sample_orf.pdf"
)
