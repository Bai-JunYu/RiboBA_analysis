#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

build_fig_s2_plot <- function(roc_df_plot, pr_df_plot, pal_method) {
  base_family <- "Arial"
  base_size <- 8
  line_w <- 0.3

  common_facet_theme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = base_size),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

  p_roc <- ggplot(roc_df_plot, aes(x = fpr, y = tpr, colour = method)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ enzyme, nrow = 1) +
    scale_colour_manual(values = pal_method, name = NULL) +
    scale_x_continuous(name = "1 - specificity", limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.2), expand = expansion(mult = 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.2), expand = expansion(mult = 0)) +
    coord_equal() +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE) +
    common_facet_theme +
    theme(legend.position = "right")

  p_pr <- ggplot(pr_df_plot, aes(x = recall, y = prec, colour = method)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ enzyme, nrow = 1) +
    scale_colour_manual(values = pal_method, name = NULL) +
    scale_x_continuous(name = "Recall", limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.2), expand = expansion(mult = 0)) +
    scale_y_continuous(name = "Precision", limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.2), expand = expansion(mult = 0)) +
    coord_equal() +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE) +
    common_facet_theme +
    theme(legend.position = "right")

  (p_roc / p_pr) +
    plot_layout(heights = c(1, 1), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "right")
}

run_fig_s2 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig_s2(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig_s2_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S2_cds_simulation_.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  req <- c("roc_df_plot", "pr_df_plot", "pal_method")
  miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))

  p <- build_fig_s2_plot(env$roc_df_plot, env$pr_df_plot, env$pal_method)
  ggsave(output_pdf, p, width = 120, height = 90, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) run_fig_s2()
