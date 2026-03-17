#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(ggplotify)
})

build_fig_s6_plot <- function(p_bias, p_corr, plig, base_family = "Arial", base_size = 9) {
  p5 <- plig$p5 +
    theme_nar(base_size = base_size, base_family = base_family, legend_inside = FALSE) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0, face = "plain", size = base_size),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )

  p3 <- plig$p3 +
    theme_nar(base_size = base_size, base_family = base_family, legend_inside = FALSE) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0, face = "plain", size = base_size),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )

  row1 <- (p_bias + plot_spacer()) + plot_layout(ncol = 2, widths = c(1, 3))
  g_corr <- ggplotify::as.grob(p_corr)
  p_corr0 <- cowplot::ggdraw() +
    cowplot::draw_grob(g_corr, x = -0.22, y = 0, width = 1.12, height = 1, clip = "off") +
    theme_void()

  blockC <- wrap_elements(full = (p5 / p3) + plot_layout(heights = c(1, 1)))

  (row1 / p_corr0 / blockC) +
    plot_layout(heights = c(0.5, 1, 1.35)) +
    plot_annotation(tag_levels = "A")
}

run_fig_s6 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig_s6(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig_s6_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S6_bias_corr.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  req <- c("p_bias", "p_corr", "plig")
  miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))

  p <- build_fig_s6_plot(env$p_bias, env$p_corr, env$plig)
  ggsave(output_pdf, p, width = 178, height = 200, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) run_fig_s6()
