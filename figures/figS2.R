#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS2_core.R",
  run_fun = "run_fig_s2",
  required = c("roc_df_plot", "pr_df_plot", "pal_method"),
  default_rds = "figure_ready_data/figS2_for_plot.rds",
  output_pdf = "figures/Fig_S2_cds_simulation_.pdf"
)
