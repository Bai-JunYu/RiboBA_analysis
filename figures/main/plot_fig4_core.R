#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(scales)
})

build_fig4_plot <- function(p_scatter2, p, p_leg) {
  shared_legend <- cowplot::get_legend(p_leg)

  p_scatter2_noleg <- p_scatter2 +
    theme(legend.position = "none", plot.margin = margin(t = 2, r = 2, b = 10, l = 8, unit = "pt"))
  p_noleg <- p +
    theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 2, l = 8, unit = "pt"))

  pA <- cowplot::ggdraw(p_scatter2_noleg) +
    cowplot::draw_label("A", x = 0, y = 1, hjust = 0, vjust = 1, size = 14, fontfamily = "Arial")
  pB <- cowplot::ggdraw(p_noleg) +
    cowplot::draw_label("B", x = 0, y = 1, hjust = 0, vjust = 1, size = 14, fontfamily = "Arial")

  core_panels <- cowplot::plot_grid(pA, pB, ncol = 1, rel_heights = c(1, 0.8), align = "v", axis = "l")
  cowplot::plot_grid(core_panels, shared_legend, ncol = 1, rel_heights = c(1, 0.10))
}

build_fig4_from_data <- function(dfp, df_mean) {
  tool_lv <- c("RiboBA", "RiboTISH", "ORF-RATER", "RiboCode", "PRICE", "RibORF")
  bio_lv <- c("uORF", "uoORF", "ncRNA", "intORF", "dORF", "doORF")
  types_keep <- c("uORF", "uoORF", "ncRNA", "intORF")
  tool_col <- c(
    "RiboBA" = "#0072B2",
    "RiboCode" = "#009E73",
    "RibORF" = "#E69F00",
    "RiboTISH" = "#D55E00",
    "PRICE" = "#CC79A7",
    "ORF-RATER" = "#C8B100"
  )[tool_lv]

  dfp <- dfp %>%
    mutate(
      tool = factor(as.character(tool), levels = tool_lv),
      biotype = factor(as.character(biotype), levels = bio_lv)
    )

  p_scatter2 <- ggplot(dfp, aes(x = x_log, y = jac_sim, colour = tool)) +
    geom_point(size = 1.4, alpha = 0.95) +
    facet_grid(pair ~ biotype, scales = "free") +
    scale_colour_manual(values = tool_col, name = NULL) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 2),
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      expand = expansion(mult = c(0.05, 0.05)),
      name = "Shared nORFs"
    ) +
    scale_y_continuous(
      name = "Jaccard similarity",
      breaks = seq(0, 0.6, 0.2),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "plain"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      axis.title.x = element_text(margin = margin(t = 2, unit = "pt")),
      axis.title.y = element_text(margin = margin(r = 2, unit = "pt"))
    ) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(shape = 15, size = 4, alpha = 1, linewidth = 0)))

  df_mean <- df_mean %>% mutate(
    tool = factor(as.character(tool), levels = tool_lv),
    orf_biotype = factor(as.character(orf_biotype), levels = types_keep)
  )

  p <- ggplot(df_mean, aes(x = idx, y = ms_mean, colour = tool, group = tool)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~orf_biotype, ncol = 4, scales = "free_x") +
    scale_colour_manual(values = tool_col, drop = FALSE, name = NULL) +
    scale_y_continuous(
      name = "MS-validated fraction (%)",
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 0.06),
      breaks = seq(0, 0.06, 0.02),
      expand = expansion(mult = c(0.01, 0.06))
    ) +
    scale_x_continuous(
      limits = c(0, NA),
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      name = "Top-ranked nORFs",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    theme_nar(base_size = 8, base_family = "Arial", line_w = 0.3, legend_inside = FALSE) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.key.width = unit(9, "pt"),
      legend.key.height = unit(8, "pt"),
      axis.title.y = element_text(margin = margin(r = 5, unit = "pt")),
      axis.text.y = element_text(margin = margin(r = 2, unit = "pt"))
    )

  leg_df <- data.frame(tool = factor(tool_lv, levels = tool_lv), x = 1, y = 1)
  p_leg <- ggplot(leg_df, aes(x, y, fill = tool)) +
    geom_tile() +
    scale_fill_manual(values = tool_col, drop = FALSE, name = NULL) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.key.width = unit(9, "pt"),
      legend.key.height = unit(8, "pt")
    )

  build_fig4_plot(p_scatter2, p, p_leg)
}

run_fig4 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig4(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig4_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_4_MS.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)

  if (exists("dfp", envir = env, inherits = FALSE) && exists("df_mean", envir = env, inherits = FALSE)) {
    out <- build_fig4_from_data(env$dfp, env$df_mean)
  } else {
    req <- c("p_scatter2", "p", "p_leg")
    miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
    if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))
    out <- build_fig4_plot(env$p_scatter2, env$p, env$p_leg)
  }

  ggsave(output_pdf, out, width = 178, height = 120, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(out)
}

if (sys.nframe() == 0) run_fig4()
