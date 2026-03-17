#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS9_core.R",
  run_fun = "run_fig_s9",
  required = c("cnt", "cnt_sc2", "sample_cols"),
  default_rds = "figure_ready_data/figS9_for_plot.rds",
  output_pdf = "figures/Fig_S9_ORF_summary.pdf"
)
