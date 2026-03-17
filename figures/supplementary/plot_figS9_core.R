#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

build_fig_s9_plot <- function(cnt, cnt_sc2, sample_cols) {
  p <- ggplot(cnt, aes(x = orf_biotype, y = N, fill = fill_grp)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c(Highlighted = "#D55E00", Other = "grey40"), guide = "none") +
    facet_wrap(~sample, ncol = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL, y = "Number of ORFs") +
    theme_nar(legend_inside = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  p_sc <- ggplot(cnt_sc2, aes(x = frac, y = start_codon, fill = sample_show)) +
    geom_col(position = position_dodge2(width = 0.85, preserve = "single"), width = 0.8, colour = NA) +
    scale_fill_manual(values = sample_cols, drop = FALSE) +
    labs(x = "Frequency", y = NULL, fill = NULL) +
    theme_nar(legend_inside = FALSE) +
    theme(axis.text.y = element_text(size = 7), plot.margin = grid::unit(c(3, 3, 3, 3), "mm"))

  design <- "
AAA
B##
"
  (p + p_sc) +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "A")
}

run_fig_s9 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) project_dir <- readline("Please input ORF_calling directory: ")
    else stop("project_dir is required. Example: run_fig_s9(project_dir='.')")
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig_s9_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S9_ORF_summary.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  req <- c("cnt", "cnt_sc2", "sample_cols")
  miss <- req[!vapply(req, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) stop("Missing objects in input RData: ", paste(miss, collapse = ", "))

  out <- build_fig_s9_plot(env$cnt, env$cnt_sc2, env$sample_cols)
  ggsave(output_pdf, out, width = 178, height = 130, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(out)
}

if (sys.nframe() == 0) run_fig_s9()
