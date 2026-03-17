#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS5_core.R",
  run_fun = "run_fig_s5",
  required = c("par_lst_split"),
  default_rds = "figure_ready_data/figS5_for_plot.rds",
  output_pdf = "figures/Fig_S5_par_split.pdf"
)
