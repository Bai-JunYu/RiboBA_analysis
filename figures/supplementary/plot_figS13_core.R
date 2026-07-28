#!/usr/bin/env Rscript

# Supplementary Figure S13: simulated codon-level P-site assignment accuracy and
# codon-level abundance reconstruction (RiboBA vs riboWaltz).
#
# Self-contained plotting from the packaged figS13_for_plot.rds. The input list
# provides:
#   codon_accuracy          : list(replicate_data, summary_data, significance_data)
#   codon_abundance_density : list(plot_data, correlation_data)

run_fig_s13 <- function(project_dir, input_rdata, output_pdf) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(grid)
  })

  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  e <- new.env()
  load(input_rdata, envir = e)
  codon_accuracy <- get("codon_accuracy", envir = e)
  codon_abundance_density <- get("codon_abundance_density", envir = e)

  condition_levels <- c("sim_rnase_i", "sim_mnase", "sim_p1")
  condition_labels <- c(
    sim_rnase_i = "RNase I",
    sim_mnase = "MNase",
    sim_p1 = "P1 nuclease"
  )
  tool_levels <- c("RiboBA", "riboWaltz")
  tool_colors <- c(RiboBA = "#0072B2", riboWaltz = "#D55E00")

  fmt_p <- function(p_value) {
    ifelse(is.na(p_value), "P = NA", ifelse(p_value < 0.001, "P < 0.001", sprintf("P = %.3f", p_value)))
  }

  tag_panel <- function(plot, tag) {
    ggdraw(plot) +
      draw_label(tag, x = 0.006, y = 0.994, hjust = 0, vjust = 1, fontface = "bold", size = 14, fontfamily = "Arial")
  }

  codon_ready <- as.data.table(copy(codon_accuracy$replicate_data))
  codon_summary <- as.data.table(copy(codon_accuracy$summary_data))
  codon_sig <- as.data.table(copy(codon_accuracy$significance_data))
  codon_ready[, tool := factor(tool, levels = tool_levels)]
  codon_summary[, tool := factor(tool, levels = tool_levels)]
  codon_ready[, condition := factor(condition, levels = condition_levels)]
  codon_summary[, condition := factor(condition, levels = condition_levels)]
  codon_ready[, condition_label := factor(condition_labels[as.character(condition)], levels = condition_labels[condition_levels])]
  codon_summary[, condition_label := factor(condition_labels[as.character(condition)], levels = condition_labels[condition_levels])]
  codon_sig[, condition := factor(condition, levels = condition_levels)]
  codon_sig[, p_label := fmt_p(p_value)]
  codon_sig[, y_position := pmin(103, y_position)]
  codon_sig[, y_bracket := y_position - 3.2]
  codon_sig[, y_tick := y_bracket - 1.2]

  p_acc <- ggplot(codon_summary, aes(x = bar_x, y = mean_codon_accuracy_percent, fill = tool)) +
    geom_col(width = 0.24, color = "grey25", linewidth = 0.25) +
    geom_errorbar(
      aes(
        ymin = mean_codon_accuracy_percent - sd_codon_accuracy_percent,
        ymax = mean_codon_accuracy_percent + sd_codon_accuracy_percent
      ),
      width = 0.08,
      linewidth = 0.32,
      na.rm = TRUE
    ) +
    geom_point(
      data = codon_ready,
      aes(x = point_x, y = codon_accuracy_percent, fill = tool),
      shape = 21,
      size = 1.45,
      stroke = 0.22,
      color = "grey15",
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = codon_sig,
      aes(x = x_min, xend = x_max, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE, linewidth = 0.3, color = "grey15"
    ) +
    geom_segment(
      data = codon_sig,
      aes(x = x_min, xend = x_min, y = y_tick, yend = y_bracket),
      inherit.aes = FALSE, linewidth = 0.3, color = "grey15"
    ) +
    geom_segment(
      data = codon_sig,
      aes(x = x_max, xend = x_max, y = y_tick, yend = y_bracket),
      inherit.aes = FALSE, linewidth = 0.3, color = "grey15"
    ) +
    geom_text(
      data = codon_sig,
      aes(x = condition_x, y = y_position, label = p_label),
      inherit.aes = FALSE,
      size = 6.3 / ggplot2::.pt,
      color = "grey10",
      family = "Arial"
    ) +
    scale_fill_manual(values = tool_colors, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(condition_levels),
      labels = unname(condition_labels[condition_levels]),
      expand = expansion(mult = c(0.07, 0.07))
    ) +
    scale_y_continuous(
      limits = c(0, 105),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(x = NULL, y = "Codon-level accuracy", fill = NULL) +
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.key.width = unit(9, "pt"),
      legend.key.height = unit(8, "pt"),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(8, 5, 4, 12)
    )

  density_dt <- as.data.table(copy(codon_abundance_density$plot_data))
  cor_dt <- as.data.table(copy(codon_abundance_density$correlation_data))
  density_dt[, condition := factor(condition, levels = condition_levels)]
  density_dt[, condition_label := factor(condition_labels[as.character(condition)], levels = condition_labels[condition_levels])]
  density_dt[, tool := factor(tool, levels = tool_levels)]
  cor_dt[, condition := factor(condition, levels = condition_levels)]
  cor_dt[, condition_label := factor(condition_labels[as.character(condition)], levels = condition_labels[condition_levels])]
  cor_dt[, tool := factor(tool, levels = tool_levels)]

  axis_min <- 0
  axis_max <- 6
  cor_dt[, label := sprintf("Spearman = %.3f", spearman_log2)]
  cor_dt[, `:=`(
    label_x = axis_min + 0.06 * (axis_max - axis_min),
    label_y = axis_max - 0.11 * (axis_max - axis_min)
  )]

  p_density <- ggplot(density_dt, aes(x = truth_log2, y = pred_log2)) +
    stat_density_2d(
      aes(
        fill = after_stat(ndensity),
        alpha = after_stat(ifelse(ndensity < 0.02, 0, sqrt(ndensity)))
      ),
      geom = "raster",
      contour = FALSE,
      n = 180,
      interpolate = TRUE,
      show.legend = FALSE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.32, color = "grey20") +
    geom_text(
      data = cor_dt,
      aes(x = label_x, y = label_y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 6.2 / ggplot2::.pt,
      color = "grey10",
      family = "Arial"
    ) +
    facet_grid(tool ~ condition_label) +
    coord_equal(xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max), expand = FALSE) +
    scale_fill_gradientn(
      colors = c("#352A87", "#2C7BB6", "#00A6CA", "#00CCBC", "#90EB9D", "#FFFF8C", "#F9D057", "#F46D43", "#D7191C"),
      na.value = "transparent"
    ) +
    scale_alpha_continuous(range = c(0, 0.95), guide = "none") +
    labs(
      x = "log2(true codon abundance + 1)",
      y = "log2(predicted codon abundance + 1)"
    ) +
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 7.6),
      axis.text = element_text(size = 6.8),
      panel.spacing = unit(2.1, "mm"),
      plot.margin = margin(2, 5, 4, 12)
    )

  out_plot <- plot_grid(
    tag_panel(p_acc, "A"),
    tag_panel(p_density, "B"),
    ncol = 2,
    rel_widths = c(1, 2.3),
    align = "h",
    axis = "tb"
  )

  ggsave(output_pdf, out_plot, width = 178, height = 90, units = "mm",
         device = grDevices::cairo_pdf, family = "Arial", dpi = 600, bg = "white")
  message("Wrote figure: ", output_pdf)
  invisible(output_pdf)
}
