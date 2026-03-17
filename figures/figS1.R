#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS1_core.R",
  run_fun = "run_fig_s1",
  required = c("sim_par_lst", "par_lst"),
  default_rds = "figure_ready_data/figS1_for_plot.rds",
  output_pdf = "figures/Fig_S1_Model_infer_simulation.pdf"
)
