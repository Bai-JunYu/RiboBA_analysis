#!/usr/bin/env Rscript

# Supplementary Figure S16: empirical P-site periodicity in human and Drosophila
# Ribo-seq libraries (RiboBA vs riboWaltz). Panel A: metagene P-site profiles
# around start/stop codons; Panel B: CDS frame-0 P-site fraction (dumbbell).
#
# Self-contained plotting from the packaged figS16_for_plot.rds. The input list
# provides:
#   human : list(meta, counts, stats)
#   fly   : list(meta, counts, stats)

run_fig_s16 <- function(project_dir, input_rdata, output_pdf) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(grid)
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
  human <- get("human", envir = e)
  fly <- get("fly", envir = e)

  tool_levels <- c("RiboBA", "riboWaltz")
  tool_colors <- c(RiboBA = "#0072B2", riboWaltz = "#D55E00")
  region_short <- c("Distance from start (nt)" = "Start codon", "Distance from stop (nt)" = "Stop codon")
  META_XLIM_START <- c(-25L,  50L)   # start codon column
  META_XLIM_STOP  <- c(-50L,  25L)   # stop codon column

  human_short_labels <- c(
    "Lucas, RNase I r1"        = "Lucas r1",
    "Lucas, RNase I r2"        = "Lucas r2",
    "Martinez, RNase I LoRes"  = "Martinez LoRes",
    "Martinez, RNase I MedRes" = "Martinez MedRes",
    "Martinez, RNase I HiRes"  = "Martinez HiRes",
    "Calviello, RNase I"       = "Calviello",
    "Darnell, MNase"           = "Darnell MNase",
    "Lucas, P1 r1"             = "Lucas P1 r1",
    "Lucas, P1 r2"             = "Lucas P1 r2"
  )
  fly_short_labels <- c(
    "Dunn JG"       = "Dunn",
    "Greenblatt EJ" = "Greenblatt",
    "Zhang H"       = "Zhang"
  )

  tag_panel <- function(plot, tag) {
    ggdraw(plot) +
      draw_label(tag, x = 0.006, y = 0.994, hjust = 0, vjust = 1,
                 fontface = "bold", size = 14, fontfamily = "Arial")
  }

  # Merge human + fly, offsetting fly plot_order so all 12 samples sort correctly.
  merge_combined <- function(hdata, fdata) {
    h_max_order <- max(hdata$meta$plot_order, na.rm = TRUE)

    dt_h <- as.data.table(copy(hdata$meta))
    dt_h[, sample_short := human_short_labels[sample_label]]
    dt_h[, plot_order_g := plot_order]
    dt_h[, species := "Human"]

    dt_f <- as.data.table(copy(fdata$meta))
    dt_f[, sample_short := fly_short_labels[sample_label]]
    dt_f[, plot_order_g := plot_order + h_max_order]
    dt_f[, species := "Drosophila"]

    st_h <- as.data.table(copy(hdata$stats))
    st_h[, sample_short := human_short_labels[sample_label]]
    st_h[, plot_order_g := plot_order]
    st_h[, species := "Human"]

    st_f <- as.data.table(copy(fdata$stats))
    st_f[, sample_short := fly_short_labels[sample_label]]
    st_f[, plot_order_g := plot_order + h_max_order]
    st_f[, species := "Drosophila"]

    list(
      meta  = rbindlist(list(dt_h, dt_f), fill = TRUE),
      stats = rbindlist(list(st_h, st_f), fill = TRUE)
    )
  }

  build_meta_plot <- function(combined, show_legend = FALSE) {
    dt <- as.data.table(copy(combined$meta))
    dt[, tool := factor(tool, levels = tool_levels)]

    sample_ord <- unique(dt[order(plot_order_g), .(sample_short)])
    dt[, sample_short := factor(sample_short, levels = sample_ord$sample_short)]
    dt[, region_label := factor(region_short[as.character(region)], levels = unname(region_short))]

    dt <- rbindlist(list(
      dt[region_label == "Start codon" & distance >= META_XLIM_START[1] & distance <= META_XLIM_START[2]],
      dt[region_label == "Stop codon"  & distance >= META_XLIM_STOP[1]  & distance <= META_XLIM_STOP[2]]
    ))

    lines3nt <- rbindlist(list(
      data.table(region_label = "Start codon", line = seq(3L,  META_XLIM_START[2], 3L)),
      data.table(region_label = "Stop codon",  line = seq(-2L, META_XLIM_STOP[1],  -3L))
    ))
    lines3nt <- lines3nt[
      (region_label == "Start codon" & line >= META_XLIM_START[1] & line <= META_XLIM_START[2]) |
      (region_label == "Stop codon"  & line >= META_XLIM_STOP[1]  & line <= META_XLIM_STOP[2])
    ]
    lines3nt[, region_label := factor(region_label, levels = unname(region_short))]

    linered <- data.table(
      region_label = factor(unname(region_short), levels = unname(region_short)),
      line = c(0, 1)
    )

    p <- ggplot(dt, aes(x = distance, y = y_percent, color = tool)) +
      geom_vline(data = lines3nt, aes(xintercept = line), inherit.aes = FALSE,
                 linetype = 3, color = "grey75", linewidth = 0.2) +
      geom_vline(data = linered, aes(xintercept = line), inherit.aes = FALSE,
                 color = "#B22222", linewidth = 0.38) +
      geom_line(linewidth = 0.5, na.rm = TRUE) +
      facet_grid(sample_short ~ region_label, scales = "free") +
      scale_x_continuous(breaks = seq(-75, 75, 25),
                         expand = expansion(mult = c(0.01, 0.01))) +
      scale_color_manual(values = tool_colors, drop = FALSE, name = NULL) +
      labs(x = "Position (nt)", y = "P-site frequency (%)", color = NULL) +
      theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
      theme(
        legend.position    = "none",
        panel.grid.major.x = element_blank(),
        axis.text.x        = element_text(size = 6.5),
        axis.text.y        = element_text(size = 6.5),
        strip.text.x       = element_text(size = 7.5),
        strip.text.y       = element_text(size = 6.5, angle = 0, hjust = 0.5,
                                          margin = margin(l = 2, r = 3)),
        plot.margin        = margin(4, 2, 4, 8)
      )

    if (show_legend) {
      p <- p + theme(
        legend.position      = "top",
        legend.direction     = "horizontal",
        legend.justification = "center",
        legend.key.width     = unit(9, "pt"),
        legend.key.height    = unit(8, "pt"),
        legend.text          = element_text(size = 7)
      )
    }
    p
  }

  build_dumbbell_plot <- function(combined) {
    stats_dt <- as.data.table(copy(combined$stats))
    setorder(stats_dt, plot_order_g)
    # rev() so first sample appears at top of y-axis
    stats_dt[, sample_short := factor(sample_short, levels = rev(unique(stats_dt$sample_short)))]

    comparison_pts <- rbindlist(list(
      stats_dt[, .(sample_short, tool = "RiboBA",    frame0 = riboba_frame0_percent)],
      stats_dt[, .(sample_short, tool = "riboWaltz", frame0 = ribowaltz_frame0_percent)]
    ))
    comparison_pts[, tool         := factor(tool, levels = tool_levels)]
    comparison_pts[, sample_short := factor(sample_short, levels = levels(stats_dt$sample_short))]

    ggplot() +
      geom_vline(xintercept = 33.333, linetype = 3, color = "grey70", linewidth = 0.28) +
      geom_segment(
        data = stats_dt,
        aes(x = ribowaltz_frame0_percent, xend = riboba_frame0_percent,
            y = sample_short, yend = sample_short),
        color = "grey55", linewidth = 0.55
      ) +
      geom_point(
        data = comparison_pts,
        aes(x = frame0, y = sample_short, color = tool),
        size = 1.8
      ) +
      scale_x_continuous(
        limits = c(25, 90),
        breaks = c(25, 50, 75),
        labels = function(x) paste0(x, "%"),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      scale_color_manual(values = tool_colors, drop = FALSE, name = NULL) +
      labs(x = "CDS frame 0 fraction", y = NULL, color = NULL) +
      theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
      theme(
        legend.position      = "top",
        legend.direction     = "horizontal",
        legend.justification = "center",
        legend.key.width     = unit(9, "pt"),
        legend.key.height    = unit(8, "pt"),
        legend.text          = element_text(size = 7),
        axis.line.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        # y-axis labels omitted: sample names shown in panel A strip text
        axis.text.y     = element_blank(),
        plot.margin     = margin(4, 6, 4, 16)
      )
  }

  combined <- merge_combined(human, fly)

  # A (lines) and B (points) each carry their own legend above the panel.
  p_meta     <- build_meta_plot(combined, show_legend = TRUE)
  p_dumbbell <- build_dumbbell_plot(combined)

  out_plot <- plot_grid(
    tag_panel(p_meta, "A"), tag_panel(p_dumbbell, "B"),
    ncol = 2, rel_widths = c(3, 1),
    align = "h", axis = "tb"
  )

  showtext_opts(dpi = 600)
  ggsave(output_pdf, out_plot, width = 178, height = 178, units = "mm",
         device = grDevices::cairo_pdf, family = "Arial", dpi = 600, bg = "white")
  message("Wrote figure: ", output_pdf)
  invisible(output_pdf)
}
