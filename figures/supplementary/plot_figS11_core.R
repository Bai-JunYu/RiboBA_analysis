#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

add_lab_title <- function(p, ttl, base_size = 8) {
  p +
    ggplot2::labs(title = ttl) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain", size = base_size, hjust = 0.5, margin = ggplot2::margin(b = 2)),
      plot.title.position = "plot"
    )
}

make_scatter <- function(res, col_map_all, title_txt = NULL, show_legend = TRUE, base_point = 0.5) {
  norf_level <- setdiff(unique(res$plot_dt$grp), "Annotated CDSs")
  if (length(norf_level) != 1) stop("grp levels wrong")
  hi_col <- unname(col_map_all[norf_level])

  p <- ggplot2::ggplot(res$plot_dt, ggplot2::aes(x = x, y = y, colour = grp)) +
    ggplot2::geom_point(alpha = 0.65, size = base_point, shape = 16, na.rm = TRUE) +
    ggplot2::geom_point(
      data = res$plot_dt[highlight == TRUE],
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21, size = 1.4, stroke = 0.35,
      fill = hi_col, colour = hi_col, alpha = 1, na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(values = col_map_all, breaks = names(col_map_all), drop = FALSE) +
    ggplot2::labs(x = "PhyloCSF per codon", y = "PhyloP per codon", colour = NULL) +
    theme_nar(legend_inside = FALSE)

  if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
  if (!is.null(title_txt) && nzchar(title_txt)) p <- add_lab_title(p, title_txt)
  p
}

build_fig_s11_plot <- function(res1, res2, res3, res1o, res2o, res3o) {
  col_map_nonov <- c("Annotated CDSs" = "grey40", "Non-overlap nORFs" = "#009E73", "Overlap nORFs" = "#D55E00")

  res1$plot_dt[, grp := fifelse(grp == "nORFs", "Non-overlap nORFs", grp)]
  res2$plot_dt[, grp := fifelse(grp == "nORFs", "Non-overlap nORFs", grp)]
  res3$plot_dt[, grp := fifelse(grp == "nORFs", "Non-overlap nORFs", grp)]
  res1o$plot_dt[, grp := fifelse(grp == "nORFs", "Overlap nORFs", grp)]
  res2o$plot_dt[, grp := fifelse(grp == "nORFs", "Overlap nORFs", grp)]
  res3o$plot_dt[, grp := fifelse(grp == "nORFs", "Overlap nORFs", grp)]

  p_non_1 <- make_scatter(res1, col_map_nonov, "Dunn JG", show_legend = TRUE, base_point = 0.5)
  p_non_2 <- make_scatter(res2, col_map_nonov, "Greenblatt EJ", show_legend = TRUE, base_point = 0.5)

  p_ov_1 <- make_scatter(res1o, col_map_nonov, "Dunn JG", show_legend = TRUE, base_point = 0.5)
  p_ov_2 <- make_scatter(res2o, col_map_nonov, "Greenblatt EJ", show_legend = TRUE, base_point = 0.5)
  p_ov_3 <- make_scatter(res3o, col_map_nonov, "Zhang H", show_legend = FALSE, base_point = 0.5)

  rowA <- (p_non_1 + p_non_2 + guide_area() + plot_spacer()) +
    plot_layout(ncol = 4, widths = c(1, 1, 0.7, 1.05), guides = "collect") &
    theme(legend.position = "right")

  rowB <- (p_ov_1 + p_ov_2 + p_ov_3 + guide_area()) +
    plot_layout(ncol = 4, widths = c(1, 1, 1, 0.6), guides = "collect") &
    theme(legend.position = "right")

  rowA_wrap <- wrap_elements(full = rowA)
  rowB_wrap <- wrap_elements(full = rowB)

  (rowA_wrap / rowB_wrap) + plot_annotation(tag_levels = "A")
}

run_fig_s11 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig_s11(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig_s11_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S11_fly_scatter_dot.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  req <- c("res1", "res2", "res3", "res1o", "res2o", "res3o")
  miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))

  p <- build_fig_s11_plot(env$res1, env$res2, env$res3, env$res1o, env$res2o, env$res3o)
  ggsave(output_pdf, p, width = 178, height = 120, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) run_fig_s11()
