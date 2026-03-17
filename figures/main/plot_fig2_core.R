#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

build_fig2_plot <- function(p_roc1, p_pr, p_bar1) {
  p_pr <- p_pr + theme(plot.margin = margin(2, 2, 6, 2, unit = "mm"))
  left_block <- wrap_elements(full = (p_roc1 / p_pr) + plot_layout(heights = c(1, 0.5)))
  left_block <- left_block & theme(plot.margin = margin(2, 0.5, 2, 2, unit = "mm"))
  p_bar1 <- p_bar1 + theme(plot.margin = margin(2, 2, 2, 0.5, unit = "mm"))

  (left_block | p_bar1) +
    plot_layout(widths = c(5, 1.2), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")
}

run_fig2 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig2(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig2_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_2_norf_simulation.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  req <- c("p_roc1", "p_pr", "p_bar1")
  miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))

  p <- build_fig2_plot(env$p_roc1, env$p_pr, env$p_bar1)
  ggsave(output_pdf, p, width = 178, height = 90, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) run_fig2()
