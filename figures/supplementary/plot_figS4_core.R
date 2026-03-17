#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

run_fig_s4 <- function(project_dir = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) {
      project_dir <- readline("Please input ORF_calling directory: ")
    } else {
      stop("project_dir is required. Example: run_fig_s4(project_dir='.')")
    }
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S4_runtime.pdf")

  base_family <- "Arial"
  base_size <- 8
  line_w <- 0.3

  df <- data.frame(
    time_sec = c(
      82 + 5378, 3196 + 462, 94 + 686,
      474, 454 + 148 + 148 + 21, 5650,
      5378 + 212, 3196 + 462, 10088 + 6,
      5378, 3196 + 462, 2743 + 24061,
      2854 + 1204 + 1236, 75133, 358 + 27387,
      10 + 5378, 3196 + 462, 5045 + 139 + 46514 + 30 + 4689 + 1687 + 280
    ),
    tool = rep(c("PRICE", "RiboBA", "RiboCode", "RiboTISH", "RibORF", "ORF-RATER"), each = 3),
    stage = rep(c("Genome preparation", "Sample preprocessing", "ORF inference"), times = 6)
  )

  order_tools <- c("PRICE", "RiboBA", "RiboCode", "RiboTISH", "RibORF", "ORF-RATER")
  order_stages <- c("Genome preparation", "Sample preprocessing", "ORF inference")

  COL_STAGE <- c(
    "Genome preparation" = "#0072B2",
    "Sample preprocessing" = "#E69F00",
    "ORF inference" = "#009E73"
  )

  df_plot <- df %>%
    mutate(
      tool = factor(tool, levels = order_tools),
      stage = factor(stage, levels = order_stages),
      time_hr = time_sec / 3600
    )

  p <- ggplot(df_plot, aes(x = tool, y = time_hr, fill = stage)) +
    geom_col(width = 0.70, color = "black", linewidth = line_w) +
    geom_text(aes(label = sprintf("%.1f", time_hr)), vjust = -0.35, size = 2.6, family = base_family) +
    scale_fill_manual(values = COL_STAGE, guide = "none") +
    facet_wrap(~ stage, nrow = 1, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = NULL, y = "Runtime (hours)") +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "plain", size = base_size),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )

  ggsave(output_pdf, p, width = 178, height = 60, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) {
  run_fig_s4()
}
