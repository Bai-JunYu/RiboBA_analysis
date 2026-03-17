#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

plot_ecdf_cds_overlap_nonoverlap2 <- function(dt,
                                              ckdt = NULL,
                                              val_col = "PhyloCSF_mean",
                                              cds_types = c("CDS", "Ext", "Trunc"),
                                              norf_types = c("doORF", "dORF", "intORF", "ncRNA", "uoORF", "uORF"),
                                              overlap_types = c("doORF", "intORF", "uoORF"),
                                              nonoverlap_types = c("dORF", "ncRNA", "uORF"),
                                              top_n = 200,
                                              inframe_cds = TRUE,
                                              inframe_norf = FALSE) {
  dt <- data.table::as.data.table(dt)
  stopifnot(all(c("orf_biotype", "is_inframe", "rpf_num", val_col) %in% names(dt)))

  val_lab <- if (val_col == "PhyloCSF_mean") "PhyloCSF per codon" else if (val_col == "phyloP_mean") "PhyloP per codon" else val_col

  cds_lab <- "Annotated CDSs"
  ov_lab <- "Overlap nORFs"
  nov_lab <- "Non-overlap nORFs"

  ann_cds <- dt[
    orf_biotype %in% cds_types & is_inframe == inframe_cds & !is.na(get(val_col)),
    .(value = get(val_col), category = cds_lab)
  ]

  n_orfs <- dt[
    orf_biotype %in% norf_types & is_inframe == inframe_norf & !is.na(rpf_num) & !is.na(get(val_col)),
    .(orf_biotype, rpf_num, value = get(val_col))
  ][order(orf_biotype, -rpf_num)][, head(.SD, top_n), by = orf_biotype]

  n_orfs[, category := data.table::fifelse(
    orf_biotype %in% overlap_types, ov_lab,
    data.table::fifelse(orf_biotype %in% nonoverlap_types, nov_lab, NA_character_)
  )]
  n_orfs <- n_orfs[!is.na(category), .(value, category)]

  real_dt <- data.table::rbindlist(list(ann_cds, n_orfs), use.names = TRUE, fill = TRUE)
  real_dt[, category := factor(category, levels = c(cds_lab, ov_lab, nov_lab))]

  p <- ggplot2::ggplot(real_dt, ggplot2::aes(x = value, colour = category)) +
    ggplot2::stat_ecdf(linewidth = 0.3) +
    ggplot2::labs(x = val_lab, y = "Cumulative fraction of ORFs", colour = NULL)

  if (!is.null(ckdt)) {
    ckdt <- data.table::as.data.table(ckdt)
    stopifnot(all(c("orf_biotype", val_col) %in% names(ckdt)))

    ck_n <- ckdt[orf_biotype %in% norf_types & !is.na(get(val_col)), .(orf_biotype, value = get(val_col))]
    ck_n[, category := data.table::fifelse(
      orf_biotype %in% overlap_types, ov_lab,
      data.table::fifelse(orf_biotype %in% nonoverlap_types, nov_lab, NA_character_)
    )]
    ck_n <- ck_n[!is.na(category), .(value, category)]
    ck_n[, category := factor(category, levels = levels(real_dt$category))]

    ck_ecdf <- ck_n[, {
      o <- order(value)
      data.table::data.table(x = value[o], y = seq_along(o) / .N)
    }, by = category]

    p <- p + ggplot2::geom_step(
      data = ck_ecdf,
      ggplot2::aes(x = x, y = y, colour = category),
      linetype = "dashed", linewidth = 0.3, show.legend = FALSE
    )
  }

  p
}

add_lab_title <- function(p, ttl, base_size = 8) {
  p +
    ggplot2::labs(title = ttl) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain", size = base_size, hjust = 0.5, margin = ggplot2::margin(b = 2)),
      plot.title.position = "plot"
    )
}

build_fig_s10_plot <- function(orf_info_annot, ckorf_info_annot) {
  ecdf_cols <- c("Annotated CDSs" = "grey40", "Overlap nORFs" = "#D55E00", "Non-overlap nORFs" = "#009E73")

  p2 <- plot_ecdf_cds_overlap_nonoverlap2(orf_info_annot$SRR942881, ckorf_info_annot$SRR942881, "PhyloCSF_mean", top_n = 50) +
    theme_nar(legend_inside = FALSE) +
    scale_colour_manual(values = ecdf_cols, drop = FALSE)
  p2 <- add_lab_title(p2, "Dunn JG")

  p2_2 <- plot_ecdf_cds_overlap_nonoverlap2(orf_info_annot$SRR10174369, ckorf_info_annot$SRR10174369, "PhyloCSF_mean", top_n = 50) +
    theme_nar(legend_inside = FALSE) +
    scale_colour_manual(values = ecdf_cols, drop = FALSE)
  p2_2 <- add_lab_title(p2_2, "Greenblatt EJ")

  p3 <- plot_ecdf_cds_overlap_nonoverlap2(orf_info_annot$SRR942881, NULL, "phyloP_mean", top_n = 50) +
    theme_nar(legend_inside = FALSE) +
    scale_colour_manual(values = ecdf_cols, drop = FALSE) +
    theme(legend.position = "none")
  p3 <- add_lab_title(p3, "Dunn JG")

  p3_2 <- plot_ecdf_cds_overlap_nonoverlap2(orf_info_annot$SRR10174369, NULL, "phyloP_mean", top_n = 50) +
    theme_nar(legend_inside = FALSE) +
    scale_colour_manual(values = ecdf_cols, drop = FALSE) +
    theme(legend.position = "none")
  p3_2 <- add_lab_title(p3_2, "Greenblatt EJ")

  p3_3 <- plot_ecdf_cds_overlap_nonoverlap2(orf_info_annot$SRR3031124, NULL, "phyloP_mean", top_n = 50) +
    theme_nar(legend_inside = FALSE) +
    scale_colour_manual(values = ecdf_cols, drop = FALSE) +
    theme(legend.position = "none")
  p3_3 <- add_lab_title(p3_3, "Zhang H")

  rowA <- (p2 + p2_2 + guide_area()) +
    plot_layout(ncol = 3, widths = c(1, 1, 0.6), guides = "collect") &
    theme(legend.position = "right")

  rowB <- (p3 + p3_2 + p3_3) + plot_layout(ncol = 3)

  rowA_wrap <- wrap_elements(full = rowA)
  rowB_wrap <- wrap_elements(full = rowB)

  (rowA_wrap / rowB_wrap) + plot_annotation(tag_levels = "A")
}

run_fig_s10 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig_s10(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig_s10_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S10_fly_Cumulative.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  req <- c("orf_info_annot", "ckorf_info_annot")
  miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))

  p <- build_fig_s10_plot(env$orf_info_annot, env$ckorf_info_annot)
  ggsave(output_pdf, p, width = 178, height = 120, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) run_fig_s10()
