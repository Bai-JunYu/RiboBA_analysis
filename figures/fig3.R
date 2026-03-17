#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/main/plot_fig3_core.R",
  run_fun = "run_fig3",
  required = c("par_lst2", "par_lst3"),
  default_rds = "figure_ready_data/fig3_for_plot.rds",
  output_pdf = "figures/Fig_3_par_summary.pdf"
)
