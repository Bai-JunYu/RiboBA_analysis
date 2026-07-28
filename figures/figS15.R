#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS15_core.R",
  run_fun = "run_fig_s15",
  required = c("posterior", "truth", "annotation"),
  default_rds = "figure_ready_data/figS15_for_plot.rds",
  output_pdf = "figures/Fig_S15_per_read_offset_uncertainty.pdf"
)
