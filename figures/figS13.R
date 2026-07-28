#!/usr/bin/env Rscript
source("figures/helpers/rds_source_wrapper.R")
run_from_rds_source(
  script_rel = "figures/supplementary/plot_figS13_core.R",
  run_fun = "run_fig_s13",
  required = c("codon_accuracy", "codon_abundance_density"),
  default_rds = "figure_ready_data/figS13_for_plot.rds",
  output_pdf = "figures/Fig_S13_simulated_codon_level_validation.pdf"
)
