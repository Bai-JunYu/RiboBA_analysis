#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidytext)
  library(tidyr)
})

build_fig_s8_plot <- function(df_ms_cnt) {
  tool_lv <- c("RiboBA", "RiboTISH", "ORF-RATER", "RiboCode", "PRICE", "RibORF")
  biotype_lv <- c("uORF", "uoORF", "ncRNA", "intORF", "dORF", "doORF")

  df_ms_cnt2 <- df_ms_cnt %>%
    mutate(
      tool = factor(as.character(tool), levels = tool_lv),
      orf_biotype = factor(as.character(orf_biotype), levels = biotype_lv)
    ) %>%
    tidyr::complete(tool, orf_biotype, fill = list(n = 0L)) %>%
    mutate(tool_ord = tidytext::reorder_within(tool, -n, orf_biotype, fun = max))

  ggplot(df_ms_cnt2, aes(x = tool_ord, y = n)) +
    geom_col(width = 0.72, fill = "grey45") +
    facet_wrap(~orf_biotype, nrow = 2, scales = "free") +
    tidytext::scale_x_reordered() +
    labs(x = NULL, y = "MS-supported nORFs") +
    theme_nar(base_size = 8, base_family = "Arial", legend_inside = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "plain"))
}

run_fig_s8 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) {
      project_dir <- readline("Please input ORF_calling directory: ")
    } else {
      stop("project_dir is required. Example: run_fig_s8(project_dir='.')")
    }
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig_s8_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S8_MS_number.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)

  if (exists("df_ms_cnt", envir = env, inherits = FALSE)) {
    p <- build_fig_s8_plot(env$df_ms_cnt)
  } else if (exists("supp_ms_orf_keep", envir = env, inherits = FALSE)) {
    tool_lv <- c("RiboBA", "RiboTISH", "ORF-RATER", "RiboCode", "PRICE", "RibORF")
    biotype_lv <- c("uORF", "uoORF", "ncRNA", "intORF", "dORF", "doORF")
    df_ms_cnt <- env$supp_ms_orf_keep %>%
      filter(orf_biotype %in% biotype_lv) %>%
      count(tool, orf_biotype, name = "n") %>%
      mutate(tool = factor(tool, levels = tool_lv), orf_biotype = factor(orf_biotype, levels = biotype_lv))
    p <- build_fig_s8_plot(df_ms_cnt)
  } else {
    stop("Input RData must contain df_ms_cnt or supp_ms_orf_keep.")
  }

  ggsave(output_pdf, p, width = 140, height = 80, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) {
  run_fig_s8()
}
