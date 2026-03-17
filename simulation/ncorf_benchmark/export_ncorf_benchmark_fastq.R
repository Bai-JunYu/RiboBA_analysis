#!/usr/bin/env Rscript

source(file.path("figures", "helpers", "prepare_utils.R"))

run_prepare_fig3_data <- function(project_dir = NULL) {
  project_dir <- get_project_dir(project_dir)
  src <- file.path(project_dir, "figures", "data", "par_lst.RData")
  dst <- file.path(project_dir, "figures", "data", "fig3_input.RData")
  e <- load_env(src)
  save_from_env(e, c("par_lst2", "par_lst3"), dst)
}

if (sys.nframe() == 0) {
  args <- parse_prepare_args()
  run_prepare_fig3_data(args$project_dir)
}
