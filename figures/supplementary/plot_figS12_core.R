#!/usr/bin/env Rscript

# Supplementary Figure S12: simulated P-site assignment periodicity and
# exact-nucleotide offset accuracy (RiboBA vs riboWaltz).
#
# Self-contained plotting from the packaged figS12_for_plot.rds. The input list
# provides:
#   meta_profile          : metagene P-site profiles around start/stop codons
#   offset_exact_accuracy : list(replicate_data, summary_data, significance_data)

run_fig_s12 <- function(project_dir, input_rdata, output_pdf) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(grid)
  })

  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  # Objects packaged in the RDS list, loaded into this frame by the wrapper.
  e <- new.env()
  load(input_rdata, envir = e)
  meta_profile <- get("meta_profile", envir = e)
  offset_exact_accuracy <- get("offset_exact_accuracy", envir = e)

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

  meta_dt <- as.data.table(copy(meta_profile))
  meta_dt[, condition := factor(condition, levels = condition_levels)]
  meta_dt[, condition_label := factor(condition_labels[as.character(condition)], levels = condition_labels[condition_levels])]
  meta_dt[, tool := factor(tool, levels = tool_levels)]
  meta_dt[, region_label := fifelse(region == "Distance from start (nt)", "Start codon", "Stop codon")]
  meta_dt[, region_label := factor(region_label, levels = c("Start codon", "Stop codon"))]

  lines3nt <- data.table(
    region_label = rep(
      c("Start codon", "Stop codon"),
      times = c(length(seq(3, 50, 3)), length(seq(-2, -50, -3)))
    ),
    line = c(seq(3, 50, 3), rev(seq(-2, -50, -3)))
  )
  lines3nt[, region_label := factor(region_label, levels = c("Start codon", "Stop codon"))]

  p_meta <- ggplot(meta_dt, aes(x = distance, y = y_percent, color = tool, fill = tool)) +
    geom_vline(data = lines3nt, aes(xintercept = line), inherit.aes = FALSE, linetype = 3, color = "grey72", linewidth = 0.22) +
    geom_ribbon(aes(ymin = pmax(0, y_percent - y_se_percent), ymax = y_percent + y_se_percent), alpha = 0.12, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.55) +
    facet_grid(condition_label ~ region_label, scales = "free_x") +
    scale_color_manual(values = tool_colors, drop = FALSE) +
    scale_fill_manual(values = tool_colors, drop = FALSE) +
    labs(x = NULL, y = "P-site count (% of all P-sites in meta-window)", color = NULL, fill = NULL) +
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.key.width = unit(9, "pt"),
      legend.key.height = unit(8, "pt"),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 6.8),
      axis.text.y = element_text(size = 6.8),
      strip.text = element_text(size = 7.6),
      plot.margin = margin(8, 5, 4, 12)
    )

  offset_ready <- as.data.table(copy(offset_exact_accuracy$replicate_data))
  offset_summary <- as.data.table(copy(offset_exact_accuracy$summary_data))
  offset_sig <- as.data.table(copy(offset_exact_accuracy$significance_data))
  offset_ready[, tool := factor(tool, levels = tool_levels)]
  offset_summary[, tool := factor(tool, levels = tool_levels)]
  offset_ready[, condition := factor(condition, levels = condition_levels)]
  offset_summary[, condition := factor(condition, levels = condition_levels)]
  offset_ready[, condition_label := factor(condition_labels[as.character(condition)], levels = condition_labels[condition_levels])]
  offset_summary[, condition_label := factor(condition_labels[as.character(condition)], levels = condition_labels[condition_levels])]
  offset_sig[, condition := factor(condition, levels = condition_levels)]
  offset_sig[, p_label := fmt_p(p_value)]

  p_offset <- ggplot(offset_summary, aes(x = bar_x, y = mean_offset_exact_accuracy_percent, fill = tool)) +
    geom_col(width = 0.24, color = "grey25", linewidth = 0.25) +
    geom_errorbar(
      aes(
        ymin = mean_offset_exact_accuracy_percent - sd_offset_exact_accuracy_percent,
        ymax = mean_offset_exact_accuracy_percent + sd_offset_exact_accuracy_percent
      ),
      width = 0.08,
      linewidth = 0.32,
      na.rm = TRUE
    ) +
    geom_point(
      data = offset_ready,
      aes(x = point_x, y = offset_exact_accuracy_percent, fill = tool),
      shape = 21,
      size = 1.45,
      stroke = 0.22,
      color = "grey15",
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = offset_sig,
      aes(x = x_min, xend = x_max, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.3,
      color = "grey15"
    ) +
    geom_segment(
      data = offset_sig,
      aes(x = x_min, xend = x_min, y = y_tick, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.3,
      color = "grey15"
    ) +
    geom_segment(
      data = offset_sig,
      aes(x = x_max, xend = x_max, y = y_tick, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.3,
      color = "grey15"
    ) +
    geom_text(
      data = offset_sig,
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
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(x = NULL, y = "Exact nucleotide-level offset accuracy", fill = NULL) +
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.margin = margin(8, 5, 4, 12)
    )

  out_plot <- plot_grid(
    tag_panel(p_meta, "A"),
    tag_panel(p_offset, "B"),
    ncol = 2,
    rel_widths = c(2.3, 1),
    align = "h",
    axis = "tb"
  )

  ggsave(output_pdf, out_plot, width = 178, height = 100, units = "mm",
         device = grDevices::cairo_pdf, family = "Arial", dpi = 600, bg = "white")
  message("Wrote figure: ", output_pdf)
  invisible(output_pdf)
}
