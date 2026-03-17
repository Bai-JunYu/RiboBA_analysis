#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS10_core.R",
  run_fun = "run_fig_s10",
  required = c("orf_info_annot", "ckorf_info_annot"),
  default_rds = "figure_ready_data/figS10_for_plot.rds",
  output_pdf = "figures/Fig_S10_fly_Cumulative.pdf"
)
