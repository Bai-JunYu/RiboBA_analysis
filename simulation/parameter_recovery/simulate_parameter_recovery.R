#!/usr/bin/env Rscript

source(file.path("figures", "helpers", "prepare_utils.R"))

run_prepare_fig_s1_data <- function(project_dir = NULL) {
  project_dir <- get_project_dir(project_dir)
  src <- file.path(project_dir, "figures", "data", "sim_for_plot.RData")
  dst <- file.path(project_dir, "figures", "data", "fig_s1_input.RData")
  e <- load_env(src)
  save_from_env(e, c("sim_par_lst", "par_lst"), dst)
}

if (sys.nframe() == 0) {
  args <- parse_prepare_args()
  run_prepare_fig_s1_data(args$project_dir)
}
