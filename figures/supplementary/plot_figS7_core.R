#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidytext)
  library(grid)
})

build_fig_s7_plot <- function(p_can2, p_norf2) {
  (p_can2 | p_norf2) +
    plot_layout(widths = c(1, 1), guides = "keep") +
    plot_annotation(tag_levels = "A") &
    theme(plot.margin = margin(1, 2, 1, 2, "mm"))
}

build_fig_s7_from_counts <- function(df_cnt) {
  tool_lv <- c("RiboBA", "RiboTISH", "ORF-RATER", "RiboCode", "PRICE", "RibORF")
  bio_norf <- c("dORF", "doORF", "uORF", "uoORF", "ncRNA", "intORF")
  bio_can <- c("Ext", "Trunc", "CDS")

  col_norf <- c(
    "uORF" = "#0072B2",
    "uoORF" = "#56B4E9",
    "ncRNA" = "#009E73",
    "intORF" = "#E69F00",
    "dORF" = "#D55E00",
    "doORF" = "#CC79A7"
  )
  col_can <- c("Ext" = "grey68", "Trunc" = "grey55", "CDS" = "grey40")

  df_norf <- df_cnt %>%
    filter(orf_biotype %in% bio_norf) %>%
    mutate(tool = factor(tool, levels = tool_lv), orf_biotype = factor(orf_biotype, levels = bio_norf))

  df_norf2 <- df_norf %>%
    group_by(sample_lab, tool) %>%
    mutate(tool_sum = sum(n, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(tool_ord = tidytext::reorder_within(tool, tool_sum, sample_lab, fun = max))

  p_norf2 <- ggplot(df_norf2, aes(x = n, y = tool_ord, fill = orf_biotype)) +
    geom_col(width = 0.72, colour = NA) +
    facet_wrap(~sample_lab, ncol = 1, scales = "free_y") +
    tidytext::scale_y_reordered() +
    scale_fill_manual(values = col_norf, breaks = rev(bio_norf), drop = FALSE, name = NULL) +
    scale_x_continuous(limits = c(0, max(df_norf2$tool_sum, na.rm = TRUE) * 1.02), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Predicted nORFs", y = NULL) +
    theme_nar(base_size = 7, base_family = "Arial", legend_inside = FALSE) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "left",
      legend.key.width = unit(7, "pt"),
      legend.key.height = unit(7, "pt"),
      legend.spacing.x = unit(2, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  df_can <- df_cnt %>%
    filter(orf_biotype %in% bio_can) %>%
    mutate(tool = factor(tool, levels = tool_lv), orf_biotype = factor(orf_biotype, levels = bio_can))

  df_can2 <- df_can %>%
    group_by(sample_lab, tool) %>%
    mutate(tool_sum = sum(n, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(tool_ord = tidytext::reorder_within(tool, tool_sum, sample_lab, fun = max))

  p_can2 <- ggplot(df_can2, aes(x = n, y = tool_ord, fill = orf_biotype)) +
    geom_col(width = 0.72, colour = NA) +
    facet_wrap(~sample_lab, ncol = 1, scales = "free_y") +
    tidytext::scale_y_reordered() +
    scale_fill_manual(values = col_can, breaks = rev(bio_can), drop = FALSE, name = NULL) +
    scale_x_continuous(limits = c(0, max(df_can2$tool_sum, na.rm = TRUE) * 1.02), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Predicted canonical ORFs", y = NULL) +
    theme_nar(base_size = 7, base_family = "Arial", legend_inside = FALSE) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "left",
      legend.key.width = unit(7, "pt"),
      legend.key.height = unit(7, "pt"),
      legend.spacing.x = unit(2, "pt")
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  build_fig_s7_plot(p_can2, p_norf2)
}

run_fig_s7 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig_s7(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig_s7_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S7_9sample_orf.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)

  if (exists("df_cnt", envir = env, inherits = FALSE)) {
    p <- build_fig_s7_from_counts(env$df_cnt)
  } else {
    req <- c("p_can2", "p_norf2")
    miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
    if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))
    p <- build_fig_s7_plot(env$p_can2, env$p_norf2)
  }

  ggsave(output_pdf, p, width = 178, height = 200, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) run_fig_s7()
