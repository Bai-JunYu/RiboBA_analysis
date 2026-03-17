#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(precrec)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source(file.path("figures", "helpers", "prepare_utils.R"))
source(file.path("figures", "helpers", "set_theme.R"))

load_res_lst_from_source <- function(src) {
  e <- load_env(src)
  x <- e$res_lst
  if (is.environment(x) && exists("res_lst", envir = x, inherits = FALSE)) {
    get("res_lst", envir = x)
  } else if (is.list(x)) {
    x
  } else {
    stop("Cannot parse res_lst from source: ", src)
  }
}

build_roc_pr_bar_objects <- function(res_lst) {
  curves_all <- imap_dfr(res_lst, function(tool_list, tool_name) {
    imap_dfr(tool_list, function(mm, sim_name) {
      as.data.frame(mm) %>%
        filter(type %in% c("ROC", "PRC")) %>%
        transmute(sim = sim_name, tool = tool_name, curvetype = type, x = x, y = y)
    })
  })

  sim_levels <- c("sim_par_i_high", "sim_par_p1", "sim_par_mnase", "sim_par_i_add5", "sim_par_i_lig")
  sim_labels <- c(
    sim_par_i_high = "RNase I",
    sim_par_p1 = "P1",
    sim_par_mnase = "MNase",
    sim_par_i_add5 = "5' add bias",
    sim_par_i_lig = "Ligation bias"
  )

  curves_all <- curves_all %>%
    mutate(
      sim = factor(sim, levels = sim_levels, labels = sim_labels[sim_levels]),
      tool = factor(tool, levels = c("RiboBA", "RiboCode", "ORF-RATER", "RibORF", "RiboTISH", "PRICE"))
    )

  roc_df_plot <- curves_all %>%
    filter(curvetype == "ROC") %>%
    rename(fpr = x, tpr = y) %>%
    mutate(tool = factor(tool, levels = c("ORF-RATER", "RiboCode", "RibORF", "RiboTISH", "PRICE", "RiboBA")))

  pr_df_plot <- curves_all %>%
    filter(curvetype == "PRC") %>%
    rename(recall = x, prec = y) %>%
    mutate(tool = factor(tool, levels = c("ORF-RATER", "RiboCode", "RibORF", "RiboTISH", "PRICE", "RiboBA")))

  cnt_long <- imap_dfr(res_lst, function(tool_list, tool_name) {
    mat <- sapply(tool_list, function(con) {
      di <- attributes(con)$data_info
      c(di[1, 3], di[1, 4])
    })

    df <- as.data.frame(mat, check.names = FALSE)
    df$label_type <- rownames(df)

    df %>%
      pivot_longer(cols = -label_type, names_to = "sim", values_to = "n") %>%
      mutate(tool = tool_name, label_type = recode(label_type, `1` = "Negative", `2` = "Positive"))
  }) %>%
    mutate(
      tool = factor(tool, levels = c("RiboBA", "RiboTISH", "ORF-RATER", "RiboCode", "PRICE", "RibORF")),
      sim = factor(sim,
                   levels = c("sim_par_i_add5", "sim_par_i_lig", "sim_par_i_high", "sim_par_p1", "sim_par_mnase"),
                   labels = c("5' add bias", "Ligation bias", "RNase I", "P1", "MNase")),
      label_type = factor(label_type, levels = c("Positive", "Negative"))
    )

  base_family <- "Arial"
  base_size <- 8
  line_w <- 0.3

  pal_method <- c(
    "RiboBA" = "#0072B2",
    "RiboCode" = "#009E73",
    "RibORF" = "#E69F00",
    "RiboTISH" = "#D55E00",
    "PRICE" = "#CC79A7",
    "ORF-RATER" = "#C8B100"
  )

  theme_rocpr_nar <- theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "plain", size = base_size), legend.position = "right", plot.margin = margin(2, 2, 2, 2, unit = "mm"))

  roc_keep <- c("RNase I", "P1", "MNase")

  p_roc1 <- ggplot(roc_df_plot[roc_df_plot$sim %in% roc_keep, ], aes(x = fpr, y = tpr, colour = tool)) +
    geom_line(linewidth = 0.5) + facet_wrap(~sim, nrow = 1) +
    scale_colour_manual(values = pal_method, name = NULL) +
    scale_x_continuous(name = "1 - specificity", limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    coord_equal() + theme_rocpr_nar + theme(legend.position = "none")

  p_roc2 <- ggplot(roc_df_plot[!(roc_df_plot$sim %in% roc_keep), ], aes(x = fpr, y = tpr, colour = tool)) +
    geom_line(linewidth = 0.5) + facet_wrap(~sim, nrow = 1) +
    scale_colour_manual(values = pal_method, name = NULL) +
    scale_x_continuous(name = "1 - specificity", limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    coord_equal() + theme_rocpr_nar + theme(legend.position = "none")

  tool_keep <- c("RiboCode", "RiboTISH", "PRICE", "RiboBA")

  p_pr <- ggplot(pr_df_plot[(pr_df_plot$tool %in% tool_keep) & (pr_df_plot$sim %in% roc_keep), ], aes(x = recall, y = prec, colour = tool)) +
    geom_line(linewidth = 0.5) + facet_wrap(~sim, nrow = 1) +
    scale_colour_manual(values = pal_method, name = NULL) +
    scale_x_continuous(name = "Recall", limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    scale_y_continuous(name = "Precision", limits = c(0.5, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    coord_equal() + theme_rocpr_nar + theme(legend.position = "none", strip.text = element_blank(), strip.background = element_blank())

  p_pr2 <- ggplot(pr_df_plot[(pr_df_plot$tool %in% tool_keep) & !(pr_df_plot$sim %in% roc_keep), ], aes(x = recall, y = prec, colour = tool)) +
    geom_line(linewidth = 0.5) + facet_wrap(~sim, nrow = 1) +
    scale_colour_manual(values = pal_method, name = NULL) +
    scale_x_continuous(name = "Recall", limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    scale_y_continuous(name = "Precision", limits = c(0.5, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    coord_equal() + theme_rocpr_nar + theme(legend.position = "none", strip.text = element_blank(), strip.background = element_blank())

  theme_bar_nar <- theme_nar(base_size = base_size, base_family = base_family, line_w = line_w, legend_inside = FALSE) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "plain", size = base_size), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(2, 2, 2, 2, unit = "mm"))

  keep_sim <- c("RNase I", "P1", "MNase")
  cnt_long$label_type <- factor(cnt_long$label_type, levels = c("Negative", "Positive"))
  pal_tool <- pal_method
  alpha_posneg <- c("Negative" = 0.25, "Positive" = 1.0)

  p_bar1 <- ggplot(cnt_long[cnt_long$sim %in% keep_sim, ], aes(x = tool, y = n, fill = tool, alpha = label_type, group = label_type)) +
    geom_col(width = 0.65, color = "black", linewidth = line_w, position = "stack") +
    facet_wrap(~sim, nrow = 3, strip.position = "right") +
    scale_fill_manual(values = pal_tool, name = NULL) +
    scale_alpha_manual(values = alpha_posneg, name = NULL) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, order = 1, override.aes = list(alpha = 1)), alpha = guide_legend(nrow = 1, byrow = TRUE, order = 2)) +
    labs(x = NULL, y = NULL) + theme_bar_nar + theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
    scale_y_continuous(breaks = scales::breaks_width(200), expand = expansion(mult = c(0, 0.05)))

  p_bar2 <- ggplot(cnt_long[!(cnt_long$sim %in% keep_sim), ], aes(x = tool, y = n, fill = tool, alpha = label_type, group = label_type)) +
    geom_col(width = 0.65, color = "black", linewidth = line_w, position = "stack") +
    facet_wrap(~sim, nrow = 3, strip.position = "right") +
    scale_fill_manual(values = pal_tool, name = NULL) +
    scale_alpha_manual(values = alpha_posneg, name = NULL) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, order = 1, override.aes = list(alpha = 1)), alpha = guide_legend(nrow = 1, byrow = TRUE, order = 2)) +
    labs(x = NULL, y = NULL) + theme_bar_nar + theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
    scale_y_continuous(breaks = scales::breaks_width(200), expand = expansion(mult = c(0, 0.05)))

  list(p_roc1 = p_roc1, p_pr = p_pr, p_bar1 = p_bar1, p_roc2 = p_roc2, p_pr2 = p_pr2, p_bar2 = p_bar2)
}

run_prepare_fig2_data <- function(project_dir = NULL, workspace_rdata = NULL) {
  project_dir <- get_project_dir(project_dir)
  public_inputs_dir <- resolve_public_inputs_dir(project_dir)
  src <- file.path(public_inputs_dir, "data", "sim_compare", "sim_info_lst3.Rdata")
  src <- normalizePath(src, mustWork = TRUE)
  res_lst <- load_res_lst_from_source(src)
  objs <- build_roc_pr_bar_objects(res_lst)
  p_roc1 <- objs$p_roc1; p_pr <- objs$p_pr; p_bar1 <- objs$p_bar1
  out <- file.path(project_dir, "figures", "data", "fig2_input.RData")
  save(p_roc1, p_pr, p_bar1, file = out)
  message("Saved: ", normalizePath(out, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  args <- parse_prepare_args()
  run_prepare_fig2_data(args$project_dir, args$workspace_rdata)
}
