#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS3_core.R",
  run_fun = "run_fig_s3",
  required = c("p_roc2", "p_pr2", "p_bar2"),
  default_rds = "figure_ready_data/figS3_for_plot.rds",
  output_pdf = "figures/Fig_S3_norf_simulation.pdf"
)
