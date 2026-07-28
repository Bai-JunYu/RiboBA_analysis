#!/usr/bin/env Rscript

# Supplementary Figure S19: RiboBA user-facing output examples (Darnell MNase,
# SRR7073124). 4-row layout, panels A-G in reading order:
#   Row 1: A (read length)   B (cleavage bias)   C (5' addition)
#   Row 2: D (steric hindrance, full width)
#   Row 3: E (ligation bias, full width)
#   Row 4: F (biotype)   G (start+stop metagene, shared y-axis)
#
# Self-contained plotting from the packaged figS19_for_plot.rds. All model-derived
# quantities were pre-extracted, so no RiboBA model object or helper.R is needed.
# The input list provides:
#   panel_A_readlen        : data.frame(x, y)  read-length relative abundance
#   panel_B_hindrance      : data.frame(dist, h_norm, enzyme, end, end_label, ...)
#   panel_C_cutbias        : data.frame(base, bias)  cleavage base bias
#   panel_D_add5           : data.frame(addition, prob)  5' non-templated addition
#   panel_E_ligation       : list(eff_f5, eff_f3)  3-mer ligation biases
#   panel_F_biotype        : data.frame(biotype, count, pct, grp) ORF biotype counts
#   panel_G_metagene       : data.table  RiboBA P-site metagene (start/stop)

run_fig_s19 <- function(project_dir, input_rdata, output_pdf) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(cowplot)
    library(showtext)
    library(stringr)
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

  BS <- 8
  theme_nar <- function(base_size = BS) {
    theme_bw(base_size = base_size, base_family = "Arial") %+replace%
      theme(
        text               = element_text(family = "Arial"),
        panel.grid.major   = element_line(linewidth = 0.25, colour = "grey88"),
        panel.grid.minor   = element_blank(),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "plain", size = base_size),
        panel.spacing      = grid::unit(3, "mm"),
        panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        axis.line          = element_blank(),
        axis.ticks         = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks.length  = grid::unit(1.1, "mm"),
        axis.title         = element_text(size = base_size),
        axis.text          = element_text(size = base_size - 1),
        legend.title       = element_text(size = base_size),
        legend.text        = element_text(size = base_size - 1),
        plot.title         = element_blank(),
        plot.subtitle      = element_blank()
      )
  }

  tag_panel <- function(plot, tag) {
    ggdraw(plot) +
      draw_label(tag, x = 0.006, y = 0.994, hjust = 0, vjust = 1,
                 fontface = "bold", size = 12, fontfamily = "Arial")
  }

  cols_base <- c(A = "#4E79A7", C = "#59A14F", G = "#E15759", T = "#F28E2B",
                 NoAdd = "#9C755F")

  e <- new.env()
  load(input_rdata, envir = e)
  panel_A_readlen <- get("panel_A_readlen", envir = e)
  panel_B_hindrance <- as.data.frame(get("panel_B_hindrance", envir = e))
  panel_C_cutbias <- as.data.frame(get("panel_C_cutbias", envir = e))
  panel_D_add5 <- get("panel_D_add5", envir = e)
  panel_E_ligation <- get("panel_E_ligation", envir = e)
  panel_F_biotype <- get("panel_F_biotype", envir = e)
  panel_G_metagene <- as.data.table(copy(get("panel_G_metagene", envir = e)))

  # ── Panel A: read length ─────────────────────────────────────────────────────
  p_A_rlen <- ggplot(panel_A_readlen, aes(x = x, y = y)) +
    geom_line(color = "#4E79A7", linewidth = 0.5) +
    geom_point(color = "#4E79A7", size = 1.0) +
    labs(x = "Read length (nt)", y = "Relative abundance") +
    theme_nar()

  # ── Panel B: steric hindrance (rebuilt from data frame) ──────────────────────
  cols_enzyme <- c("RNase I" = "#0072B2", "MNase" = "#009E73", "P1" = "#D55E00", "Other" = "#666666")
  df_h <- panel_B_hindrance
  df_h$end_label <- factor(ifelse(df_h$end == "5prime", "5'", "3'"), levels = c("5'", "3'"))
  p_B_hin <- ggplot(
      df_h,
      aes(x = dist, y = h_norm, color = enzyme, linetype = end_label,
          group = interaction(sample, end))
    ) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_color_manual(values = cols_enzyme, name = "Enzyme") +
    scale_linetype_manual(values = c("solid", "solid")) +
    labs(x = "Distance to P-site (nt)", y = "Steric hindrance") +
    theme_nar() +
    theme(legend.position = "none")

  # ── Panel C: cleavage base bias (rebuilt from data frame) ────────────────────
  cols_cbase <- c(A = "#4E79A7", C = "#59A14F", G = "#E15759", T = "#F28E2B")
  df_cb <- panel_C_cutbias
  df_cb$base <- factor(df_cb$base, levels = c("A", "C", "G", "T"))
  p_C_cut <- ggplot(df_cb, aes(x = base, y = bias, fill = base)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = cols_cbase, guide = "none") +
    labs(x = NULL, y = "Cleavage bias") +
    theme_nar()

  # ── Panel D: 5' non-templated addition ───────────────────────────────────────
  nms_add5 <- as.character(panel_D_add5$addition)
  p_D_add5 <- ggplot(panel_D_add5, aes(x = addition, y = prob, fill = addition)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = cols_base[nms_add5], guide = "none") +
    labs(x = NULL, y = "Probability") +
    theme_nar()

  # ── Panel E: combined 5'+3' ligation bias (all 64 3-mers) ────────────────────
  eff_f5 <- setNames(panel_E_ligation$eff_f5$effect, panel_E_ligation$eff_f5$kmer)
  eff_f3 <- setNames(panel_E_ligation$eff_f3$effect, panel_E_ligation$eff_f3$kmer)
  df_lig <- rbind(
    data.frame(kmer = names(eff_f5), effect = as.numeric(eff_f5), end = "5′", stringsAsFactors = FALSE),
    data.frame(kmer = names(eff_f3), effect = as.numeric(eff_f3), end = "3′", stringsAsFactors = FALSE)
  )
  df_lig$kmer <- factor(df_lig$kmer, levels = sort(unique(df_lig$kmer)))
  df_lig$end  <- factor(df_lig$end, levels = c("5′", "3′"))
  end_cols <- c("5′" = "#D55E00", "3′" = "#0072B2")
  p_E_lig <- ggplot(df_lig, aes(x = kmer, y = effect, color = end)) +
    geom_point(size = 2.0, position = position_dodge(width = 0.4)) +
    scale_color_manual(values = end_cols, name = "Terminal") +
    labs(x = NULL, y = "Ligation bias") +
    theme_nar() +
    theme(
      axis.text.x     = element_text(size = BS - 2, family = "mono",
                                     angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "right",
      legend.key.size = grid::unit(5, "pt"),
      legend.text     = element_text(size = BS - 1),
      legend.title    = element_text(size = BS - 1)
    )

  # ── Panel F: ncORF biotype composition ───────────────────────────────────────
  df_bio <- as.data.frame(panel_F_biotype)
  df_bio_vis <- df_bio[df_bio$count > 0, ]
  # biotype is already an ordered factor (levels ascending by count, so CDS -- the
  # largest -- is the last level and appears at the top of the y-axis); keep it.
  df_bio_vis$biotype <- droplevels(df_bio_vis$biotype)
  bio_cols <- c(annotated = "#4E79A7", novel = "#F28E2B")
  p_F <- ggplot(df_bio_vis, aes(x = count, y = biotype, fill = grp)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = formatC(count, format = "d", big.mark = ",")),
              hjust = -0.15, size = 2.0, family = "Arial") +
    scale_fill_manual(values = bio_cols, guide = "none") +
    scale_x_log10(
      breaks = c(1, 10, 100, 1000, 10000),
      labels = function(x) ifelse(x >= 1000, paste0(x / 1000, "k"), as.character(x)),
      expand = expansion(mult = c(0, 0.40))
    ) +
    labs(x = "Translated ORF count (log scale)", y = NULL) +
    theme_nar() +
    theme(plot.margin = margin(4, 8, 4, 4))

  # ── Panel G: P-site metagene (start + stop, shared y) ────────────────────────
  dt_cd <- panel_G_metagene
  y_lim <- c(0, max(dt_cd$count, na.rm = TRUE) * 1.06)
  lines3nt <- rbind(
    data.table(region = "Distance from start (nt)", line = seq(3L, 50L, 3L)),
    data.table(region = "Distance from stop (nt)",  line = rev(seq(-2L, -50L, -3L)))
  )
  frame_cols <- c("0" = "#D55E00", "1" = "#009E73", "2" = "#0072B2")
  frame_offset <- list(
    "Distance from start (nt)" = 1L,
    "Distance from stop (nt)"  = 2L
  )
  frame_of <- function(distance, reg) {
    b <- frame_offset[[reg]]
    factor(as.character(((distance - b) %% 3L)), levels = c("0", "1", "2"))
  }
  make_meta <- function(reg, xlab) {
    dp <- dt_cd[region == reg]
    dp[, frame := frame_of(distance, reg)]
    l3 <- lines3nt[region == reg]
    ggplot(dp, aes(x = distance, y = count, fill = frame)) +
      geom_vline(data = l3, aes(xintercept = line), inherit.aes = FALSE,
                 linetype = 3, color = "grey55", linewidth = 0.28) +
      geom_col(width = 0.8, linewidth = 0) +
      scale_fill_manual(values = frame_cols, name = "Frame", drop = FALSE) +
      coord_cartesian(ylim = y_lim) +
      scale_y_continuous(
        labels = function(x) ifelse(x >= 1000, paste0(x / 1000, "k"), as.character(x))
      ) +
      labs(x = xlab, y = "P-site count") +
      theme_nar() +
      theme(plot.margin = margin(3, 6, 3, 4))
  }
  p_H <- make_meta("Distance from start (nt)", "Distance from start codon (nt)")
  p_I <- make_meta("Distance from stop (nt)",  "Distance from stop codon (nt)")

  frame_legend <- cowplot::get_legend(
    p_H + theme(
      legend.position = "right",
      legend.key.size = grid::unit(5, "pt"),
      legend.text     = element_text(size = BS - 1),
      legend.title    = element_text(size = BS - 1)
    )
  )
  p_H_nol <- p_H + theme(legend.position = "none")
  p_I_noy <- p_I + theme(
    legend.position = "none",
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin  = margin(3, 6, 3, 1)
  )
  p_GH_core <- plot_grid(p_H_nol, p_I_noy, nrow = 1, align = "h", axis = "tb",
                         rel_widths = c(1.2, 1))
  p_GH <- plot_grid(p_GH_core, frame_legend, nrow = 1, rel_widths = c(1, 0.14))

  # ── Final assembly ───────────────────────────────────────────────────────────
  row1 <- plot_grid(
    tag_panel(p_A_rlen, "A"), tag_panel(p_C_cut, "B"), tag_panel(p_D_add5, "C"),
    nrow = 1, align = "hv", axis = "tblr"
  )
  row2 <- tag_panel(p_B_hin, "D")
  row3 <- tag_panel(p_E_lig, "E")
  row4 <- plot_grid(
    tag_panel(p_F, "F"), tag_panel(p_GH, "G"),
    nrow = 1, rel_widths = c(1, 2), align = "hv", axis = "tblr"
  )
  fig <- plot_grid(row1, row2, row3, row4, ncol = 1, rel_heights = c(1, 0.96, 0.87, 1.0))

  showtext_opts(dpi = 600)
  ggsave(output_pdf, fig, width = 7, height = 6.5,
         device = grDevices::cairo_pdf, family = "Arial")
  message("Wrote figure: ", output_pdf)
  invisible(output_pdf)
}
