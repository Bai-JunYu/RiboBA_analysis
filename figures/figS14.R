#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS14_core.R",
  run_fun = "run_fig_s14",
  required = c("rnase_p1", "mnase"),
  default_rds = "figure_ready_data/figS14_for_plot.rds",
  output_pdf = "figures/Fig_S14_stratified_offset_probability_distributions.pdf"
)
