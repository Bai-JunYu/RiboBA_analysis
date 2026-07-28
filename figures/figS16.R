#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS16_core.R",
  run_fun = "run_fig_s16",
  required = c("human", "fly"),
  default_rds = "figure_ready_data/figS16_for_plot.rds",
  output_pdf = "figures/Fig_S16_empirical_periodicity.pdf"
)
