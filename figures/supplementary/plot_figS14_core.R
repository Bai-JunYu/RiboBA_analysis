#!/usr/bin/env Rscript

# Supplementary Figure S14: stratified offset probability distributions showing
# how protocol-induced biases shift the offset within and across read-length
# classes (ground truth vs RiboBA posterior vs riboWaltz fixed offset).
#
# Self-contained plotting from the packaged figS14_for_plot.rds. The input list
# provides:
#   rnase_p1 : list(distribution, traditional_offsets)  -- RNase I / P1 by length x 5'-addition
#   mnase    : list(distribution, traditional_offsets)  -- MNase by length x P-site-offset group

run_fig_s14 <- function(project_dir, input_rdata, output_pdf) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggridges)
    library(cowplot)
    library(grid)
  })

  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  e <- new.env()
  load(input_rdata, envir = e)
  rnase_p1 <- get("rnase_p1", envir = e)
  mnase <- get("mnase", envir = e)

  condition_levels <- c("sim_rnase_i", "sim_mnase", "sim_p1")
  condition_labels <- c(
    sim_rnase_i = "RNase I",
    sim_mnase = "MNase",
    sim_p1 = "P1 nuclease"
  )
  source_colors <- c(
    "Ground truth" = "grey72",
    "RiboBA posterior" = "#0072B2",
    "riboWaltz fixed offset" = "#D55E00"
  )
  grp_levels <- c("Low MAP offset", "Middle MAP offset", "High MAP offset")

  offset_theme <- function() {
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.key.width = unit(9, "pt"),
        legend.key.height = unit(8, "pt"),
        legend.text = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 8),
        plot.margin = margin(5, 5, 4, 8)
      )
  }

  build_rnase_p1_plot <- function(plot_dt, vlines, show_legend = FALSE) {
    truth_dt <- as.data.table(copy(plot_dt[source == "Ground truth"]))
    riboba_dt <- as.data.table(copy(plot_dt[source == "RiboBA posterior"]))

    len_order <- plot_dt[order(qwidth), unique(qwidth_label)]
    truth_dt[, qwidth_label := factor(qwidth_label, levels = len_order)]
    riboba_dt[, qwidth_label := factor(qwidth_label, levels = len_order)]

    key_cols <- c(
      "condition", "condition_label", "qwidth", "qwidth_label",
      "length_rank", "add_status", "offset"
    )
    riboba_dt <- merge(
      truth_dt[, ..key_cols],
      riboba_dt[, c(key_cols, "probability"), with = FALSE],
      by = key_cols,
      all.x = TRUE
    )
    riboba_dt[is.na(probability), probability := 0]
    riboba_dt[, qwidth_label := factor(qwidth_label, levels = len_order)]

    vdt <- as.data.table(copy(vlines[!is.na(traditional_offset)]))
    vdt[, qwidth_label := factor(qwidth_label, levels = len_order)]
    vdt[, y_num := as.integer(factor(qwidth, levels = sort(unique(qwidth)))), by = condition_label]

    ridge_scale <- 0.82
    divider_dt <- unique(truth_dt[, .(condition_label, add_status, qwidth)])[
      , .(y_div = seq_len(.N)), by = .(condition_label, add_status)
    ]

    p <- ggplot() +
      geom_segment(
        data = divider_dt,
        aes(x = -Inf, xend = Inf, y = y_div, yend = y_div),
        color = "grey80", linewidth = 0.25, inherit.aes = FALSE
      ) +
      geom_ridgeline(
        data = truth_dt,
        aes(x = offset, y = qwidth_label, height = probability, fill = "Ground truth"),
        color = "grey45", linewidth = 0.2, alpha = 0.85, scale = ridge_scale
      ) +
      geom_ridgeline(
        data = riboba_dt,
        aes(x = offset, y = qwidth_label, height = probability, color = "RiboBA posterior"),
        fill = NA, linewidth = 0.6, scale = ridge_scale, key_glyph = draw_key_path
      ) +
      geom_segment(
        data = vdt,
        aes(x = traditional_offset, xend = traditional_offset,
            y = y_num, yend = y_num + ridge_scale * 0.45,
            color = "riboWaltz fixed offset"),
        linewidth = 0.85
      ) +
      facet_grid(condition_label ~ add_status, scales = "free") +
      scale_fill_manual(values = source_colors["Ground truth"], breaks = "Ground truth", name = NULL) +
      scale_color_manual(
        values = source_colors[c("RiboBA posterior", "riboWaltz fixed offset")],
        breaks = c("RiboBA posterior", "riboWaltz fixed offset"),
        labels = c("RiboBA posterior offset", "riboWaltz fixed offset"),
        name = NULL
      ) +
      scale_y_discrete(expand = expansion(mult = c(0.02, 0.14))) +
      labs(x = "Offset from 5' end (nt)", y = NULL) +
      guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 2, override.aes = list(fill = NA))
      ) +
      offset_theme()

    if (!show_legend) p <- p + theme(legend.position = "none")
    p
  }

  build_mnase_plot <- function(plot_dt, vlines, show_legend = FALSE) {
    truth_dt <- as.data.table(copy(plot_dt[source == "Ground truth"]))
    riboba_dt <- as.data.table(copy(plot_dt[source == "RiboBA posterior"]))

    truth_dt[, map_offset_group := factor(map_offset_group, levels = grp_levels)]
    riboba_dt[, map_offset_group := factor(map_offset_group, levels = grp_levels)]

    len_order <- truth_dt[order(qwidth), unique(qwidth_label)]
    truth_dt[, qwidth_label := factor(qwidth_label, levels = len_order)]
    riboba_dt[, qwidth_label := factor(qwidth_label, levels = len_order)]

    key_cols <- c(
      "condition", "condition_label", "qwidth", "qwidth_label",
      "length_rank", "map_offset_group", "offset"
    )
    riboba_dt <- merge(
      truth_dt[, ..key_cols],
      riboba_dt[, c(key_cols, "probability"), with = FALSE],
      by = key_cols,
      all.x = TRUE
    )
    riboba_dt[is.na(probability), probability := 0]
    riboba_dt[, map_offset_group := factor(map_offset_group, levels = grp_levels)]
    riboba_dt[, qwidth_label := factor(qwidth_label, levels = len_order)]

    vdt <- as.data.table(copy(vlines[!is.na(traditional_offset) & map_offset_group %in% grp_levels]))
    vdt[, map_offset_group := factor(map_offset_group, levels = grp_levels)]
    vdt[, y_num := as.integer(factor(qwidth, levels = sort(unique(qwidth)))), by = condition_label]

    ridge_scale <- 0.82
    divider_dt <- unique(truth_dt[, .(condition_label, map_offset_group, qwidth)])[
      , .(y_div = seq_len(.N)), by = .(condition_label, map_offset_group)
    ]

    p <- ggplot() +
      geom_segment(
        data = divider_dt,
        aes(x = -Inf, xend = Inf, y = y_div, yend = y_div),
        color = "grey80", linewidth = 0.25, inherit.aes = FALSE
      ) +
      geom_ridgeline(
        data = truth_dt,
        aes(x = offset, y = qwidth_label, height = probability, fill = "Ground truth"),
        color = "grey45", linewidth = 0.2, alpha = 0.85, scale = ridge_scale
      ) +
      geom_ridgeline(
        data = riboba_dt,
        aes(x = offset, y = qwidth_label, height = probability, color = "RiboBA posterior"),
        fill = NA, linewidth = 0.6, scale = ridge_scale, key_glyph = draw_key_path
      ) +
      geom_segment(
        data = vdt,
        aes(x = traditional_offset, xend = traditional_offset,
            y = y_num, yend = y_num + ridge_scale * 0.45,
            color = "riboWaltz fixed offset"),
        linewidth = 0.85
      ) +
      facet_grid(
        condition_label ~ map_offset_group, scales = "free_x",
        labeller = labeller(map_offset_group = c(
          "Low MAP offset"    = "Short P-site offset",
          "Middle MAP offset" = "Medium P-site offset",
          "High MAP offset"   = "Long P-site offset"
        ))
      ) +
      scale_fill_manual(values = source_colors["Ground truth"], breaks = "Ground truth", name = NULL) +
      scale_color_manual(
        values = source_colors[c("RiboBA posterior", "riboWaltz fixed offset")],
        breaks = c("RiboBA posterior", "riboWaltz fixed offset"),
        labels = c("RiboBA posterior offset", "riboWaltz fixed offset"),
        name = NULL
      ) +
      scale_y_discrete(expand = expansion(mult = c(0.02, 0.14))) +
      labs(x = "Offset from 5' end (nt)", y = NULL) +
      guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 2, override.aes = list(fill = NA))
      ) +
      offset_theme()

    if (!show_legend) p <- p + theme(legend.position = "none")
    p
  }

  tag_panel <- function(plot, tag) {
    ggdraw(plot) +
      draw_label(tag, x = 0.006, y = 0.994, hjust = 0, vjust = 1,
                 fontface = "bold", size = 14, fontfamily = "Arial")
  }

  rnase_p1_dt <- as.data.table(copy(rnase_p1$distribution))
  rnase_p1_vlines <- as.data.table(copy(rnase_p1$traditional_offsets))
  mnase_dt <- as.data.table(copy(mnase$distribution))
  mnase_vlines <- as.data.table(copy(mnase$traditional_offsets))

  p_rnase_p1 <- build_rnase_p1_plot(rnase_p1_dt, rnase_p1_vlines, show_legend = FALSE)
  p_mnase <- build_mnase_plot(mnase_dt, mnase_vlines, show_legend = FALSE)
  legend <- get_legend(build_rnase_p1_plot(rnase_p1_dt, rnase_p1_vlines, show_legend = TRUE))

  # A has 2 facet columns, B has 3; A is narrowed to 2/3 width so the offset-axis
  # unit length matches between A and B, leaving the remaining third blank.
  a_row <- plot_grid(
    tag_panel(p_rnase_p1, "A"),
    ggdraw(),
    ncol = 2,
    rel_widths = c(2, 1)
  )

  out_plot <- plot_grid(
    legend,
    a_row,
    tag_panel(p_mnase, "B"),
    ncol = 1,
    rel_heights = c(0.12, 1.34, 0.74)
  )

  ggsave(output_pdf, out_plot, width = 178, height = 170, units = "mm",
         device = grDevices::cairo_pdf, family = "Arial", dpi = 600, bg = "white")
  message("Wrote figure: ", output_pdf)
  invisible(output_pdf)
}
