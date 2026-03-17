#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(scales)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(paste0("^", file_arg), args)
  if (length(idx) > 0) {
    path <- sub(file_arg, "", args[idx[1]])
    if (nzchar(path) && path != "-" && file.exists(path)) {
      return(dirname(normalizePath(path, mustWork = TRUE)))
    }
  }

  if (sys.nframe() > 0) {
    for (i in rev(seq_len(sys.nframe()))) {
      fr <- sys.frame(i)
      if (exists("ofile", envir = fr, inherits = FALSE)) {
        ofile <- get("ofile", envir = fr, inherits = FALSE)
        if (is.character(ofile) && length(ofile) > 0 && nzchar(ofile[1]) && ofile[1] != "-" && file.exists(ofile[1])) {
          return(dirname(normalizePath(ofile[1], mustWork = TRUE)))
        }
      }
    }
  }

  normalizePath(getwd(), mustWork = TRUE)
}

resolve_set_theme <- function(project_dir) {
  p <- file.path(project_dir, "figures", "helpers", "set_theme.R")
  if (!file.exists(p)) stop("Cannot find set_theme.R: ", p)
  p
}

normalize_split <- function(x) {
  if (is.null(x$split1) && !is.null(x[[1]])) x$split1 <- x[[1]]
  if (is.null(x$split2) && !is.null(x[[2]])) x$split2 <- x[[2]]
  x
}

sample_to_lab <- function(samp) {
  case_when(
    str_detect(samp, regex("RNase\\s*I", ignore_case = TRUE)) ~ "RNase I",
    str_detect(samp, regex("\\bMNase\\b", ignore_case = TRUE)) ~ "MNase",
    str_detect(samp, regex("\\bP1\\b", ignore_case = TRUE)) ~ "P1",
    TRUE ~ "Other"
  )
}

cut_prob_plot <- function(par0, eps = 1e-12) {
  p5 <- as.numeric(par0$prob_hd$p5)
  p3 <- as.numeric(par0$prob_hd$p3)
  p5 <- pmax(p5, eps)
  p3 <- pmax(p3, eps)
  tibble(pmf5 = p5, pmf3 = p3)
}

get_pmf <- function(par) {
  cp <- cut_prob_plot(par)
  list(f5 = cp$pmf5, f3 = cp$pmf3)
}

get_s7 <- function(x) {
  v <- x$cut_bias$s7[c("A", "C", "G", "T")]
  as.numeric(v)
}

cor_align <- function(a, b) {
  n <- min(length(a), length(b))
  suppressWarnings(cor(a[seq_len(n)], b[seq_len(n)], method = "pearson"))
}

get_add5 <- function(x) {
  v <- x$prob_add5
  if (length(v) == 5) names(v)[5] <- "no_add"
  as.numeric(v[c("A", "C", "G", "T", "no_add")])
}

build_fig_s5_plot <- function(par_lst_split) {
  BASE4_NAR <- c("A" = "#0072B2", "C" = "#E69F00", "G" = "#009E73", "T" = "#D55E00")

  par_split <- lapply(par_lst_split, normalize_split)

  dfA <- purrr::imap_dfr(par_split, function(s2, samp) {
    p1 <- get_pmf(s2$split1)
    p2 <- get_pmf(s2$split2)
    tibble(
      sample = samp,
      lab = sample_to_lab(samp),
      `5'` = cor_align(p1$f5, p2$f5),
      `3'` = cor_align(p1$f3, p2$f3)
    )
  }) %>%
    tidyr::pivot_longer(cols = c("5'", "3'"), names_to = "end", values_to = "r")

  target <- names(par_split)[stringr::str_detect(names(par_split), regex("MNase", ignore_case = TRUE))][1]
  if (is.na(target)) stop("No MNase sample found in par_lst_split.")

  s7_1 <- get_s7(par_split[[target]]$split1)
  s7_2 <- get_s7(par_split[[target]]$split2)
  dfB <- tibble(base = factor(c("A", "C", "G", "T"), levels = c("A", "C", "G", "T")), split1 = s7_1, split2 = s7_2)

  dfC <- purrr::imap_dfr(par_split, function(s2, samp) {
    tibble(
      sample = samp,
      lab = sample_to_lab(samp),
      end5 = cor_align(s2$split1$eff_f5, s2$split2$eff_f5),
      end3 = cor_align(s2$split1$eff_f3, s2$split2$eff_f3)
    )
  }) %>%
    tidyr::pivot_longer(cols = c(end5, end3), names_to = "end", values_to = "r") %>%
    mutate(end = recode(end, end5 = "5'", end3 = "3'"))

  a5_1 <- get_add5(par_split[[target]]$split1)
  a5_2 <- get_add5(par_split[[target]]$split2)
  dfD <- tibble(base = factor(c("A", "C", "G", "T", "no_add"), levels = c("A", "C", "G", "T", "no_add")), split1 = a5_1, split2 = a5_2)

  base_family <- "Arial"
  base_size <- 8
  line_w <- 0.3

  pA <- ggplot(dfA, aes(x = end, y = r, fill = end)) +
    geom_violin(width = 0.9, alpha = 0.5, color = NA) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.9) +
    scale_fill_manual(values = c("5'" = "#a6cee3", "3'" = "#b2df8a"), guide = "none") +
    scale_y_continuous(limits = c(0.995, 1), expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = NULL, y = "Pearson's r") +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE)

  pB <- ggplot(dfB, aes(split1, split2, color = base)) +
    geom_abline(slope = 1, intercept = 0, linewidth = line_w, linetype = "dashed", color = "grey60") +
    geom_point(size = 1.8) +
    scale_color_manual(values = BASE4_NAR, name = NULL) +
    labs(x = "Cut bias (split 1)", y = "Cut bias (split 2)") +
    coord_equal(expand = TRUE) +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE)

  pC <- ggplot(dfC, aes(x = end, y = r, fill = end)) +
    geom_violin(width = 0.9, alpha = 0.5, color = NA) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.9) +
    scale_fill_manual(values = c("5'" = "#a6cee3", "3'" = "#b2df8a"), guide = "none") +
    scale_y_continuous(limits = c(0.994, 1), expand = expansion(mult = c(0.1, 0.1))) +
    labs(x = NULL, y = "Pearson's r") +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE)

  COL_ADD5 <- c(BASE4_NAR, "no_add" = "#666666")
  pD <- ggplot(dfD, aes(split1, split2, color = base)) +
    geom_abline(slope = 1, intercept = 0, linewidth = line_w, linetype = "dashed", color = "grey60") +
    geom_point(size = 1.8) +
    scale_color_manual(values = COL_ADD5, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.1, 0.1))) +
    labs(x = "Addition probability (split 1)", y = "Addition probability (split 2)") +
    coord_equal(expand = TRUE) +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE)

  row1 <- (pA | plot_spacer() | pC) + plot_layout(widths = c(0.8, 0.2, 1))
  row2 <- (pB | plot_spacer() | pD) + plot_layout(widths = c(0.8, 0.2, 1))

  (row1 / row2) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "A")
}

run_fig_s5 <- function(project_dir = NULL, input_rdata = NULL, output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) {
      project_dir <- readline("Please input ORF_calling directory: ")
    } else {
      stop("project_dir is required. Example: run_fig_s5(project_dir='.')")
    }
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  source(resolve_set_theme(project_dir))

  if (is.null(input_rdata)) input_rdata <- file.path(project_dir, "figures", "data", "fig_s5_input.RData")
  if (is.null(output_pdf)) output_pdf <- file.path(project_dir, "figures", "Fig_S5_par_split.pdf")

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  if (!exists("par_lst_split", envir = env)) stop("Input RData must contain par_lst_split.")

  p <- build_fig_s5_plot(env$par_lst_split)
  ggsave(output_pdf, p, width = 178, height = 120, units = "mm", device = grDevices::cairo_pdf, dpi = 600, bg = "white")
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) {
  run_fig_s5()
}
