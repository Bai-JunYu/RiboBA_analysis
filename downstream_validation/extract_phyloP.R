#!/usr/bin/env Rscript

source(file.path("figures", "helpers", "prepare_utils.R"))

run_prepare_fig_s10_data <- function(project_dir = NULL, workspace_rdata = NULL) {
  project_dir <- get_project_dir(project_dir)
  public_inputs_dir <- resolve_public_inputs_dir(project_dir)
  src <- file.path(public_inputs_dir, "data", "fly", "phyloCSF.RData")
  src <- normalizePath(src, mustWork = TRUE)

  e <- load_env(src)
  save_from_env(e, c("orf_info_annot", "ckorf_info_annot"), file.path(project_dir, "figures", "data", "fig_s10_input.RData"))
}

if (sys.nframe() == 0) {
  args <- parse_prepare_args()
  run_prepare_fig_s10_data(args$project_dir, args$workspace_rdata)
}
