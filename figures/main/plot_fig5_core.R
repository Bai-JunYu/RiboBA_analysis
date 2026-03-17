#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

build_fig5_plot <- function(p2_3, p_non_3, p_all) {
  p_non_3 <- p_non_3 +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key.width = grid::unit(10, "pt"),
      legend.key.height = grid::unit(6, "pt"),
      legend.spacing.x = grid::unit(2, "pt")
    )

  top_row <- (p2_3 + plot_spacer() + p_non_3) +
    plot_layout(ncol = 3, widths = c(1.3, 0.005, 0.8))

  top_row <- top_row & theme(plot.margin = grid::unit(c(2, 2, 2, 2), "mm"))

  (top_row / wrap_elements(full = p_all)) +
    plot_layout(heights = c(0.95, 1.05)) +
    plot_annotation(tag_levels = "A")
}

run_fig5 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig5(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig5_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_5_fly.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  req <- c("p2_3", "p_non_3", "p_all")
  miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))

  p <- build_fig5_plot(env$p2_3, env$p_non_3, env$p_all)
  ggsave(output_pdf, p, width = 178, height = 140, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) run_fig5()
