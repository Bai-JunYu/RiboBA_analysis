#!/usr/bin/env Rscript

source(file.path("figures", "helpers", "prepare_utils.R"))

run_prepare_fig_s2_data <- function(project_dir = NULL, workspace_rdata = NULL) {
  project_dir <- get_project_dir(project_dir)
  public_inputs_dir <- resolve_public_inputs_dir(project_dir)
  src <- file.path(public_inputs_dir, "data", "simulate1", "simu_fastq", "curves_all.RData")
  src <- normalizePath(src, mustWork = TRUE)

  e <- load_env(src)
  ce <- e$curves_all
  if (is.environment(ce)) {
    curves_df <- get("curves_all", envir = ce)
  } else {
    curves_df <- ce
  }

  stopifnot(is.data.frame(curves_df))
  req <- c("enzyme", "method", "curvetype")
  miss <- req[!req %in% names(curves_df)]
  if (length(miss) > 0) stop("Missing columns in curves_all: ", paste(miss, collapse = ", "))

  roc_df_plot <- subset(curves_df, toupper(curvetype) %in% c("ROC", "ROCS"))
  if (!"fpr" %in% names(roc_df_plot) || !"tpr" %in% names(roc_df_plot)) {
    stop("ROC data must contain fpr/tpr columns.")
  }

  pr_df_plot <- subset(curves_df, toupper(curvetype) %in% c("PR", "PRC", "PRCS"))
  if (!"recall" %in% names(pr_df_plot) || !"prec" %in% names(pr_df_plot)) {
    stop("PR data must contain recall/prec columns.")
  }

  method_levels <- c("RiboCode", "RibORF", "RiboTISH", "PRICE", "RiboBA")
  roc_df_plot$method <- factor(as.character(roc_df_plot$method), levels = method_levels)
  pr_df_plot$method <- factor(as.character(pr_df_plot$method), levels = method_levels)
  roc_df_plot <- roc_df_plot[!is.na(roc_df_plot$method), , drop = FALSE]
  pr_df_plot <- pr_df_plot[!is.na(pr_df_plot$method), , drop = FALSE]

  pal_method <- c(
    "RiboBA" = "#0072B2",
    "RiboCode" = "#009E73",
    "RibORF" = "#E69F00",
    "RiboTISH" = "#D55E00",
    "PRICE" = "#CC79A7"
  )
  pal_method <- pal_method[method_levels]

  out <- file.path(project_dir, "figures", "data", "fig_s2_input.RData")
  save(roc_df_plot, pr_df_plot, pal_method, file = out)
  message("Saved: ", normalizePath(out, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  args <- parse_prepare_args()
  run_prepare_fig_s2_data(args$project_dir, args$workspace_rdata)
}
