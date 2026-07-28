#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(scales)
})

# Pitfall 2 fix: cairo_pdf mis-kerns TrueType fonts (e.g. Arial) on this Linux
# setup even when the font is correctly installed. Load Arial via showtext so
# text is shaped/outlined correctly; keep cairo_pdf as the ggsave device.
# dpi must match the ggsave dpi (600) for correct text sizing.
if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
  .arial_dir <- "/usr/share/fonts/truetype/msttcorefonts"
  if (file.exists(file.path(.arial_dir, "Arial.ttf")) &&
      !"Arial" %in% sysfonts::font_families()) {
    sysfonts::font_add(
      family     = "Arial",
      regular    = file.path(.arial_dir, "Arial.ttf"),
      bold       = file.path(.arial_dir, "Arial_Bold.ttf"),
      italic     = file.path(.arial_dir, "Arial_Italic.ttf"),
      bolditalic = file.path(.arial_dir, "Arial_Bold_Italic.ttf")
    )
  }
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 600)
}

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

get_project_dir <- function(project_dir = NULL) {
  if (is.null(project_dir) || !nzchar(project_dir)) {
    if (interactive()) {
      project_dir <- readline("Please input ORF_calling directory: ")
    } else {
      stop(
        "project_dir is required.\n",
        "Example: run_fig_s1(project_dir = '.')"
      )
    }
  }
  if (!dir.exists(project_dir)) {
    stop("Invalid project_dir: ", project_dir)
  }
  if (!dir.exists(file.path(project_dir, "figures"))) {
    stop("Invalid ORF_calling directory (missing 'figures'): ", project_dir)
  }
  normalizePath(project_dir, mustWork = TRUE)
}

resolve_helper_plot <- function(project_dir, helper_plot = NULL) {
  candidates <- c(
    helper_plot,
    file.path(project_dir, "figures", "helpers", "helper_plot.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop(
      "Cannot find helper_plot.R. Checked:\n",
      paste0("- ", candidates, collapse = "\n")
    )
  }
  normalizePath(hit, mustWork = TRUE)
}

prepare_fig_s1_data <- function(sim_par_lst, par_lst) {
  conds_default <- c(
    "RNase I (high)", "RNase I (low)", "MNase", "P1", "5'add-base", "Ligation model"
  )
  names(sim_par_lst) <- conds_default
  names(par_lst) <- conds_default

  pmf_df_from_par <- function(par_obj, cond_label) {
    pmf_mat <- cut_prob_plot(par_obj)
    L5 <- par_obj$ribo_size[1]
    S5 <- par_obj$ribo_size[2]
    S3 <- par_obj$ribo_size[3]
    L3 <- par_obj$ribo_size[4]

    pmf5 <- as.numeric(pmf_mat[, "pmf5"])
    pmf3 <- as.numeric(pmf_mat[, "pmf3"])
    dist5 <- -(S5 + (L5:1))
    dist3 <- S3 + 0:(L3 - 1)

    bind_rows(
      tibble(cond = cond_label, end = "5'", dist = dist5, pmf = pmf5),
      tibble(cond = cond_label, end = "3'", dist = dist3, pmf = pmf3)
    ) %>% mutate(end = factor(end, levels = c("5'", "3'")))
  }

  mean_named_vec <- function(vlist) {
    keys <- sort(unique(unlist(lapply(vlist, names))))
    mat <- sapply(vlist, function(v) v[keys])
    rowMeans(mat, na.rm = TRUE)
  }

  susc_ratio <- function(bias_named) {
    v <- as.numeric(bias_named[c("A", "C", "G", "T")])
    names(v) <- c("A", "C", "G", "T")
    v <- pmin(pmax(v, 1e-8), 1 - 1e-8)
    s <- -log(v)
    s / min(s)
  }

  conds_a <- conds_default[1:4]
  panel_a_sim <- map2(sim_par_lst[conds_a], conds_a, pmf_df_from_par) %>%
    bind_rows() %>%
    mutate(cond = factor(cond, levels = conds_a))

  panel_a_est <- imap(par_lst[conds_a], function(reps, nm) {
    map(reps, ~ pmf_df_from_par(.x, cond_label = nm)) %>% bind_rows(.id = "rep")
  }) %>%
    bind_rows() %>%
    group_by(cond, end, dist) %>%
    summarise(mean = mean(pmf), sd = sd(pmf), .groups = "drop") %>%
    mutate(cond = factor(cond, levels = conds_a))

  cond6 <- conds_default[6]
  sim6 <- sim_par_lst[[cond6]]
  reps6 <- par_lst[[cond6]]

  est6_mean_f5 <- mean_named_vec(lapply(reps6, `[[`, "eff_f5"))
  est6_mean_f3 <- mean_named_vec(lapply(reps6, `[[`, "eff_f3"))

  k5 <- intersect(names(sim6$eff_f5), names(est6_mean_f5))
  k3 <- intersect(names(sim6$eff_f3), names(est6_mean_f3))

  panel_b_df <- bind_rows(
    tibble(end = "5'", sim = unname(sim6$eff_f5[k5]), est = unname(est6_mean_f5[k5])),
    tibble(end = "3'", sim = unname(sim6$eff_f3[k3]), est = unname(est6_mean_f3[k3]))
  ) %>% mutate(end = factor(end, levels = c("5'", "3'")))

  panel_b_ann <- panel_b_df %>%
    group_by(end) %>%
    summarise(
      r = suppressWarnings(cor(sim, est, method = "pearson", use = "complete.obs")),
      p = suppressWarnings(cor.test(sim, est, method = "pearson")$p.value),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("r = %.3f, p = %.2g", r, p))

  corr_per_rep <- function(reps_vec_list, sim_vec) {
    sapply(reps_vec_list, function(v) suppressWarnings(cor(sim_vec, v, method = "pearson")))
  }

  panel_c_df <- bind_rows(
    tibble(end = "5'", r = corr_per_rep(lapply(reps6, `[[`, "eff_f5"), sim6$eff_f5)),
    tibble(end = "3'", r = corr_per_rep(lapply(reps6, `[[`, "eff_f3"), sim6$eff_f3))
  ) %>% mutate(end = factor(end, levels = c("5'", "3'")))

  cond5 <- conds_default[5]
  sim5 <- sim_par_lst[[cond5]]$prob_add5
  names(sim5)[5] <- "No add"
  lev5 <- c("A", "C", "G", "T", "No add")

  est5_tbl <- map(par_lst[[cond5]], ~ .x$prob_add5) %>% do.call(rbind, .) %>% as.data.frame()
  colnames(est5_tbl) <- names(sim5)

  panel_d_sum <- as_tibble(est5_tbl) %>%
    pivot_longer(everything(), names_to = "base", values_to = "prob") %>%
    group_by(base) %>%
    summarise(mean = mean(prob), sd = sd(prob), .groups = "drop") %>%
    mutate(base = factor(base, levels = lev5))

  panel_d_sim <- tibble(
    base = factor(names(sim5), levels = lev5),
    prob = as.numeric(sim5)
  )

  cond3 <- conds_default[3]
  sim3_bias <- sim_par_lst[[cond3]]$cut_bias$s7

  panel_e_sim <- tibble(
    base = factor(names(susc_ratio(sim3_bias)), levels = c("A", "C", "G", "T")),
    value = as.numeric(susc_ratio(sim3_bias)),
    type = "Ground truth"
  )

  est3_list <- lapply(par_lst[[cond3]], function(p) susc_ratio(p$cut_bias$s7))
  panel_e_sum <- bind_rows(lapply(seq_along(est3_list), function(i) {
    tibble(
      rep = i,
      base = factor(names(est3_list[[i]]), levels = c("A", "C", "G", "T")),
      value = as.numeric(est3_list[[i]])
    )
  })) %>%
    group_by(base) %>%
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")

  list(
    conds_default = conds_default,
    panel_a_sim = panel_a_sim,
    panel_a_est = panel_a_est,
    panel_b_df = panel_b_df,
    panel_b_ann = panel_b_ann,
    panel_c_df = panel_c_df,
    panel_d_sim = panel_d_sim,
    panel_d_sum = panel_d_sum,
    panel_e_sim = panel_e_sim,
    panel_e_sum = panel_e_sum
  )
}

plot_fig_s1 <- function(fig_data, base_family = "Arial", base_size = 8, line_w = 0.3) {
  conds_a <- fig_data$conds_default[1:4]
  pal_a <- c(
    "RNase I (high)" = "#0072B2",
    "RNase I (low)" = "#E69F00",
    "MNase" = "#009E73",
    "P1" = "#D55E00"
  )

  series_cols_a <- c(
    "Ground truth (simulated)" = "#0072B2",
    "RiboBA estimated (mean ± SD)" = "black"
  )
  p_a <- ggplot() +
    geom_line(
      data = fig_data$panel_a_sim,
      aes(dist, pmf, colour = "Ground truth (simulated)"),
      linewidth = 0.4
    ) +
    geom_errorbar(
      data = fig_data$panel_a_est,
      aes(dist, ymin = mean - sd, ymax = mean + sd,
          colour = "RiboBA estimated (mean ± SD)"),
      width = 0.4,
      linewidth = 0.4
    ) +
    facet_grid(rows = vars(cond), cols = vars(end), scales = "free_x") +
    scale_colour_manual(name = NULL, values = series_cols_a) +
    scale_x_continuous(
      breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 2),
      minor_breaks = NULL,
      labels = label_number(accuracy = 1, style_negative = "minus")
    ) +
    labs(x = "Distance to P-site (nt)", y = "Cut probability") +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w) +
    theme(legend.position = "right", legend.key.width = grid::unit(12, "pt"))

  p_b <- ggplot(fig_data$panel_b_df, aes(x = sim, y = est)) +
    geom_point(alpha = 0.85, size = 1) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
    facet_wrap(~end, nrow = 1) +
    geom_text(
      data = fig_data$panel_b_ann,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.1,
      vjust = 1.1,
      size = 7 / 2.845,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
    scale_y_continuous(labels = label_number(style_negative = "minus")) +
    labs(x = "Ground-truth ligation efficiency", y = "Estimated ligation efficiency") +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w)

  fill_map <- c("5'" = nar_palette()[6], "3'" = nar_palette()[3])
  p_c <- ggplot(fig_data$panel_c_df, aes(end, r, fill = end)) +
    geom_violin(width = 0.9, alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8, linewidth = line_w) +
    geom_jitter(width = 0.08, size = 1, alpha = 0.8) +
    scale_fill_manual(values = fill_map, guide = "none") +
    labs(x = NULL, y = "Ligation efficiency correlation (r)") +
    coord_cartesian(ylim = c(0.88, 0.95)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w)

  base_cols_d <- c(
    "A" = nar_palette()[1],
    "C" = nar_palette()[8],
    "G" = nar_palette()[3],
    "T" = nar_palette()[2],
    "No add" = "grey60"
  )
  p_d <- ggplot() +
    geom_errorbar(
      data = fig_data$panel_d_sum,
      aes(base, ymin = pmax(mean - sd, 0), ymax = mean + sd,
          colour = "RiboBA estimated (mean ± SD)"),
      width = 0.8,
      linewidth = 0.4
    ) +
    geom_point(
      data = fig_data$panel_d_sim,
      aes(base, prob, colour = "Ground truth (simulated)"),
      shape = 17,
      size = 1.6
    ) +
    scale_colour_manual(name = NULL, values = series_cols_a) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = NULL, y = "5' addition probability") +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  base_cols_e <- c(
    "A" = nar_palette()[1],
    "C" = nar_palette()[8],
    "G" = nar_palette()[3],
    "T" = nar_palette()[2]
  )
  p_e <- ggplot() +
    geom_errorbar(
      data = fig_data$panel_e_sum,
      aes(base, ymin = pmax(mean - sd, 0), ymax = mean + sd,
          colour = "RiboBA estimated (mean ± SD)"),
      width = 0.8,
      linewidth = 0.4
    ) +
    geom_point(
      data = fig_data$panel_e_sim,
      aes(base, value, colour = "Ground truth (simulated)"),
      shape = 17,
      size = 1.6
    ) +
    scale_colour_manual(name = NULL, values = series_cols_a) +
    labs(x = NULL, y = "Relative cleavability") +
    theme_nar(base_size = base_size, base_family = base_family, line_w = line_w)

  row1 <- p_a
  row2 <- (p_b + p_c) + plot_layout(ncol = 2, widths = c(2.1, 0.6))
  row3 <- (p_e + p_d) + plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
    theme(legend.position = "right")

  (row1 / row2 / row3) +
    plot_layout(heights = c(2.5, 1.45, 1.0)) +
    plot_annotation(tag_levels = "A")
}

run_fig_s1 <- function(
    project_dir = NULL,
    helper_plot = NULL,
    input_rdata = NULL,
    output_pdf = NULL) {
  project_dir <- get_project_dir(project_dir)
  helper_plot <- resolve_helper_plot(project_dir, helper_plot)
  source(file.path(project_dir, "figures", "helpers", "set_theme.R"))
  if (is.null(input_rdata)) {
    input_rdata <- file.path(project_dir, "figures", "data", "fig_s1_input.RData")
  }
  if (is.null(output_pdf)) {
    output_pdf <- file.path(project_dir, "figures", "Fig_S1_Model_infer_simulation.pdf")
  }
  source(helper_plot)

  if (!file.exists(input_rdata)) {
    stop("Input data file not found: ", input_rdata)
  }

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  if (!exists("sim_par_lst", envir = env) || !exists("par_lst", envir = env)) {
    stop("RData must contain both 'sim_par_lst' and 'par_lst'.")
  }

  fig_data <- prepare_fig_s1_data(env$sim_par_lst, env$par_lst)
  p <- plot_fig_s1(fig_data)

  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = output_pdf,
    plot = p,
    width = 178,
    height = 230,
    units = "mm",
    device = grDevices::cairo_pdf,
    dpi = 600
  )

  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) {
  run_fig_s1()
}
