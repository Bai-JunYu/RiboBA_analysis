#!/usr/bin/env Rscript

source(file.path("simulation", "ncorf_benchmark", "build_donor_profiles.R"))
source(file.path("figures", "helpers", "prepare_utils.R"))

run_prepare_fig_s3_data <- function(project_dir = NULL, workspace_rdata = NULL) {
  project_dir <- get_project_dir(project_dir)
  public_inputs_dir <- resolve_public_inputs_dir(project_dir)
  src <- file.path(public_inputs_dir, "data", "sim_compare", "sim_info_lst3.Rdata")
  src <- normalizePath(src, mustWork = TRUE)
  res_lst <- load_res_lst_from_source(src)
  objs <- build_roc_pr_bar_objects(res_lst)
  p_roc2 <- objs$p_roc2; p_pr2 <- objs$p_pr2; p_bar2 <- objs$p_bar2
  out <- file.path(project_dir, "figures", "data", "fig_s3_input.RData")
  save(p_roc2, p_pr2, p_bar2, file = out)
  message("Saved: ", normalizePath(out, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  args <- parse_prepare_args()
  run_prepare_fig_s3_data(args$project_dir, args$workspace_rdata)
}
