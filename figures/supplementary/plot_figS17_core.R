#!/usr/bin/env Rscript

# Supplementary Figure S17: paired library-level comparison of the area under the
# Fig. 4B MS-validation curve between RiboBA and each competing tool, per ncORF
# biotype. Points are per-library area differences (RiboBA - competitor); the
# crossbar is the median; annotations are BH-corrected one-sided paired Wilcoxon
# q-values.
#
# Self-contained plotting from the packaged figS17_for_plot.rds. The input list
# provides:
#   deltas : per-library area differences (columns: competitor, display_biotype,
#            delta, auc_source, ...)
#   tests  : paired-test summary with bh_q_wilcoxon_greater per biotype x competitor

run_fig_s17 <- function(project_dir, input_rdata, output_pdf) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(scales)
    library(showtext)
  })

  showtext_auto()
  .afd <- "/usr/share/fonts/truetype/msttcorefonts"
  if (file.exists(file.path(.afd, "Arial.ttf")) && !"Arial" %in% sysfonts::font_families()) {
    font_add("Arial",
             regular    = file.path(.afd, "Arial.ttf"),
             bold       = file.path(.afd, "arialbd.ttf"),
             italic     = file.path(.afd, "ariali.ttf"),
             bolditalic = file.path(.afd, "arialbi.ttf"))
  }

  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  e <- new.env()
  load(input_rdata, envir = e)
  delta_plot <- as.data.table(copy(get("deltas", envir = e)))
  test_plot <- as.data.table(copy(get("tests", envir = e)))

  tool_levels <- c("RiboBA", "RiboTISH", "ORF-RATER", "RiboCode", "PRICE", "RibORF")
  competitor_levels <- setdiff(tool_levels, "RiboBA")
  display_levels <- c("uORF", "uoORF", "lncORF", "intORF")
  tool_cols <- c(
    RiboBA = "#0072B2",
    RiboTISH = "#D55E00",
    `ORF-RATER` = "#C8B100",
    RiboCode = "#009E73",
    PRICE = "#CC79A7",
    RibORF = "#E69F00"
  )

  format_q <- function(x) {
    ifelse(
      is.na(x),
      "NA",
      ifelse(x < 0.001, formatC(x, format = "e", digits = 1), formatC(x, format = "f", digits = 3))
    )
  }

  delta_plot[, `:=`(
    competitor = factor(competitor, levels = competitor_levels),
    display_biotype = factor(display_biotype, levels = display_levels),
    auc_source = factor(
      auc_source,
      levels = c("displayed_grid_auc", "full_raw_auc_fallback"),
      labels = c("Displayed grid AUC", "Full raw AUC fallback")
    )
  )]
  test_plot[, `:=`(
    competitor = factor(competitor, levels = competitor_levels),
    display_biotype = factor(display_biotype, levels = display_levels),
    figure_label = paste0("q=", format_q(bh_q_wilcoxon_greater))
  )]

  delta_range <- delta_plot[, .(
    ymin = min(delta, na.rm = TRUE),
    ymax = max(delta, na.rm = TRUE)
  ), by = display_biotype]
  delta_range[, label_y := ymax + pmax(1.5, 0.12 * (ymax - ymin))]
  label_pos <- merge(test_plot, delta_range[, .(display_biotype, label_y)], by = "display_biotype", all.x = TRUE)
  median_delta <- delta_plot[, .(median_delta = median(delta, na.rm = TRUE)), by = .(display_biotype, competitor)]

  p <- ggplot(delta_plot, aes(x = competitor, y = delta, colour = competitor)) +
    geom_hline(yintercept = 0, linewidth = 0.28, colour = "grey35") +
    geom_point(
      shape = 16,
      position = position_jitter(width = 0.13, height = 0, seed = 1),
      size = 1.35,
      alpha = 0.72
    ) +
    geom_crossbar(
      data = median_delta,
      aes(x = competitor, y = median_delta, ymin = median_delta, ymax = median_delta),
      width = 0.52,
      linewidth = 0.35,
      colour = "black",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = label_pos,
      aes(x = competitor, y = label_y, label = figure_label),
      colour = "black",
      size = 2.0,
      lineheight = 0.86,
      vjust = 0,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ display_biotype, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = tool_cols[competitor_levels], drop = FALSE, guide = "none") +
    scale_shape_manual(values = c("Displayed grid AUC" = 16, "Full raw AUC fallback" = 17),
                       drop = FALSE, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.30))) +
    labs(
      x = NULL,
      y = expression(Delta * " area under MS-validation curve (relative to competitor)")
    ) +
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      panel.spacing.x = unit(5, "mm"),
      panel.grid = element_blank()
    )

  showtext_opts(dpi = 600)
  ggsave(output_pdf, p, width = 178, height = 120, units = "mm",
         device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Wrote figure: ", output_pdf)
  invisible(output_pdf)
}
