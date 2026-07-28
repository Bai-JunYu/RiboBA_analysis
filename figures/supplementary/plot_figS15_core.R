#!/usr/bin/env Rscript

# Supplementary Figure S15: RiboBA per-read offset posterior vs ground truth for
# three nucleases. For each nuclease, three reads of the same length carry
# different ground-truth offsets (sequence-context driven); RiboBA's posterior
# tracks the true offset as it shifts, which a fixed length-specific offset
# cannot.
#
# Self-contained plotting from the packaged figS15_for_plot.rds. The input list
# provides:
#   posterior  : per-read RiboBA posterior over candidate offsets
#   truth      : ground-truth offset distribution (simulated, known)
#   annotation : per-panel read-length labels

run_fig_s15 <- function(project_dir, input_rdata, output_pdf) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(grid)
  })

  # cairo_pdf mis-kerns TrueType Arial on some machines (text metrics vs glyph
  # mismatch). Load Arial via showtext so layout metrics match the drawn glyphs.
  # dpi must equal the ggsave dpi (600).
  if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
    .afd <- "/usr/share/fonts/truetype/msttcorefonts"
    if (file.exists(file.path(.afd, "Arial.ttf")) && !"Arial" %in% sysfonts::font_families()) {
      sysfonts::font_add("Arial",
        regular = file.path(.afd, "Arial.ttf"),
        bold = file.path(.afd, "Arial_Bold.ttf"),
        italic = file.path(.afd, "Arial_Italic.ttf"),
        bolditalic = file.path(.afd, "Arial_Bold_Italic.ttf"))
    }
    showtext::showtext_auto(); showtext::showtext_opts(dpi = 600)
  }

  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  e <- new.env()
  load(input_rdata, envir = e)
  dt <- as.data.table(copy(get("posterior", envir = e)))
  truth <- as.data.table(copy(get("truth", envir = e)))
  peak_dt <- as.data.table(copy(get("annotation", envir = e)))

  enzyme_levels <- c("RNase I", "MNase", "P1 nuclease")
  dt[, enzyme_label := factor(enzyme_label, levels = enzyme_levels)]
  truth[, enzyme_label := factor(enzyme_label, levels = enzyme_levels)]

  x_max <- max(dt$offset, truth$offset)

  truth_fill <- "#BDBDBD"      # ground truth: neutral grey area
  post_fill  <- "#E6820D"      # RiboBA posterior: orange bars (visual protagonist)
  line_col   <- "grey55"       # read-intrinsic cleavage likelihood: recede to background

  p <- ggplot() +
    geom_col(
      data = truth,
      aes(x = offset, y = truth_percent, fill = "Ground-truth offset (simulated)"),
      width = 0.92, color = NA
    ) +
    geom_line(
      data = dt,
      aes(x = offset, y = base_percent, linetype = "RiboBA offset (before P-site prior)"),
      color = line_col, linewidth = 0.7
    ) +
    geom_col(
      data = dt,
      aes(x = offset, y = y_percent, fill = "RiboBA posterior offset"),
      width = 0.5, color = NA
    ) +
    geom_text(
      data = peak_dt, aes(x = x_max, y = Inf, label = lab),
      hjust = 1, vjust = 1.3, size = 6.5 / .pt, lineheight = 0.9, color = "grey30", inherit.aes = FALSE
    ) +
    facet_grid(enzyme_label ~ read_col, switch = "y") +
    scale_x_continuous(breaks = scales::breaks_width(3),
                       labels = scales::label_number(accuracy = 1, style_negative = "minus")) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(
      name = NULL,
      values = c("Ground-truth offset (simulated)" = truth_fill, "RiboBA posterior offset" = post_fill),
      breaks = c("Ground-truth offset (simulated)", "RiboBA posterior offset")
    ) +
    scale_linetype_manual(name = NULL, values = c("RiboBA offset (before P-site prior)" = "22")) +
    labs(x = "Offset from 5' end of read (nt)", y = "Offset probability (%)") +
    guides(fill = guide_legend(order = 1),
           linetype = guide_legend(order = 2, override.aes = list(color = line_col))) +
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
    theme(
      legend.position = "top", legend.direction = "horizontal", legend.justification = "center",
      legend.key.width = unit(14, "pt"), legend.key.height = unit(8, "pt"),
      legend.text = element_text(size = 6.5), legend.margin = margin(0, 0, 2, 0),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(6, "pt"), panel.spacing.y = unit(3, "pt"),
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.5, face = "plain"),
      strip.text.x = element_text(size = 6.5, face = "plain"), strip.placement = "outside",
      axis.title = element_text(size = 6.5), axis.text = element_text(size = 6.5),
      plot.margin = margin(6, 6, 4, 6)
    )

  ggsave(output_pdf, p, width = 178, height = 135, units = "mm",
         device = grDevices::cairo_pdf, family = "Arial", dpi = 600, bg = "white")
  message("Wrote figure: ", output_pdf)
  invisible(output_pdf)
}
