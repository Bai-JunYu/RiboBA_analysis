library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)

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

resolve_set_theme <- function() {
  candidates <- c(
    file.path(get_script_dir(), "helpers", "set_theme.R"),
    file.path(get_script_dir(), "figures", "helpers", "set_theme.R"),
    file.path(getwd(), "figures", "helpers", "set_theme.R"),
    file.path(dirname(get_script_dir()), "figures", "helpers", "set_theme.R")
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("Cannot find set_theme.R. Checked:\n", paste0("- ", candidates, collapse = "\n"))
  }
  normalizePath(hit, mustWork = TRUE)
}

source(resolve_set_theme())
# Fig 2 9 samples par overview ####
plot_hd_nar <- function(
    par_lst,
    p0 = 0.9,
    eps = 1e-3,
    base_family = "Arial",
    base_size = 7
) {
  lst_par <- par_lst
  
  # ---- 1) par -> tidy ----
  one_sample <- function(par0, id) {
    stopifnot(!is.null(par0$prob_hd$p5), !is.null(par0$prob_hd$p3))
    d5 <- tibble(sample = id, end = "5′", k = seq_along(par0$prob_hd$p5), m = as.numeric(par0$prob_hd$p5))
    d3 <- tibble(sample = id, end = "3′", k = seq_along(par0$prob_hd$p3), m = as.numeric(par0$prob_hd$p3))
    bind_rows(d5, d3)
  }
  
  df_h <- imap_dfr(lst_par, one_sample) %>%
    mutate(
      m = pmax(m, eps),
      h = -log(m)
    ) %>%
    group_by(sample, end) %>%
    mutate(
      h_norm = h / (-log(p0)),
      dist = if_else(end == "5′", -rev(row_number()), row_number())
    ) %>%
    ungroup() %>%
    mutate(
      end    = factor(end, levels = c("5′", "3′")),
      sample = factor(sample, levels = names(lst_par))
    )
  
  # ---- 2) ribo_size dist shift (according to your original formula) ----
  ribo_size <- sapply(lst_par, function(x) x$ribo_size)  # matrix-like
  tmp_v1 <- -(ribo_size[2, ])  # 5' shift
  tmp_v2 <-  (ribo_size[3, ])  # 3' shift
  
  is5 <- df_h$end == "5′"
  df_h$dist[is5]  <- df_h$dist[is5]  + tmp_v1[match(df_h$sample[is5], names(tmp_v1))]
  df_h$dist[!is5] <- df_h$dist[!is5] + tmp_v2[match(df_h$sample[!is5], names(tmp_v2))]
  
  # ---- 3) Color grouping (enzyme) ----
  df_plot <- df_h %>%
    mutate(
      enzyme = case_when(
        str_detect(sample, regex("RNase\\s*I", ignore_case = TRUE)) ~ "RNase I",
        str_detect(sample, regex("Martinez\\s*,\\s*(Epicentre|TruSeq)", ignore_case = TRUE)) ~ "RNase I",
        str_detect(sample, regex("\\bMNase\\b", ignore_case = TRUE)) ~ "MNase",
        str_detect(sample, regex("\\bP1\\b", ignore_case = TRUE)) ~ "P1",
        TRUE ~ NA_character_
      ),
      enzyme = factor(enzyme, levels = c("RNase I", "MNase", "P1"))
    )
  
  # Add a dist = 0 breakpoint for each sample × end (no dots, just to break the line in the middle more natural)
  df_zero <- df_plot %>%
    distinct(sample, end, enzyme) %>%
    mutate(dist = 0, h_norm = NA_real_, k = NA_integer_, m = NA_real_, h = NA_real_)
  
  df_plot2 <- bind_rows(df_plot, df_zero)
  
  cols_enzyme <- c("RNase I" = "#0072B2", "MNase" = "#009E73", "P1" = "#D55E00")
  
  ggplot(df_plot2,
         aes(x = dist, y = h_norm,
             colour = enzyme,
             group = interaction(sample, end))) +
    geom_line(linewidth = 0.4, na.rm = TRUE) +
    geom_point(size = 1, stroke = 0.15, na.rm = TRUE) +
    scale_colour_manual(values = cols_enzyme, name = NULL, drop = FALSE) +
    scale_x_continuous(
      breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 4),
      minor_breaks = NULL
    ) +
    scale_y_continuous(
      position = "right",
      breaks = seq(0, 65, by = 30),
      limits = c(0, 65)
    ) +
    labs(x = "Distance to P-site (nt)", y = "Relative steric hindrance") +
    facet_grid(
      rows = vars(sample), cols = vars(end),
      scales = "free_x", switch = "y"
    ) +
    theme_nar(base_size = base_size, base_family = base_family, legend_inside = FALSE) +
    theme(
      legend.position = "none",
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, vjust = 0.5),
      # Hindrance maps generally do not require a strong mesh (you can keep it if you want)
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88"),
      panel.grid.major.x = element_blank()
    )
}

BASE4_NAR <- c(
  "A" = "#0072B2",
  "C" = "#E69F00",
  "G" = "#009E73",
  "T" = "#D55E00"
)

plot_mnase_base_specificity_nar <- function(
    par_0,
    normalize = c("sum","mean","none"),
    base_family = "Arial",
    base_size = 7,
    legend_inside = TRUE,
    legend_pos = c(0.78, 0.98)
) {
  normalize <- match.arg(normalize)
  
  m <- par_0$cut_bias$s7[c("A","C","G","T")]
  m <- pmin(pmax(as.numeric(m), 1e-8), 1 - 1e-8)
  names(m) <- c("A","C","G","T")
  
  haz <- -log(m)
  rel <- switch(normalize,
                sum  = haz / sum(haz),
                mean = haz / mean(haz),
                none = haz / min(haz))
  
  df <- tibble(
    base  = factor(names(rel), levels = c("A","C","G","T")),
    value = as.numeric(rel)
  )
  
  ggplot(df, aes(base, value, fill = base)) +
    geom_col(width = 0.70, colour = "black", linewidth = 0.4) +
    scale_fill_manual(values = BASE4_NAR, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = NULL, y = "MNase cut bias", fill = NULL) +
    theme_nar(base_size = base_size, base_family = base_family,
              legend_inside = legend_inside, legend_pos = legend_pos) +
    theme(
      legend.position = "none",
      axis.ticks.x = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
}


plot_add5_prob_nar <- function(
    par_5add,
    normalize = TRUE,
    base_family = "Arial",
    base_size = 7,
    legend_inside = TRUE,
    legend_pos = c(0.78, 0.98)
) {
  req <- c("A","C","G","T","No add")
  if (!all(req %in% names(par_5add))) stop("par_5add must have names: A,C,G,T,No add")
  
  v <- as.numeric(par_5add[req])
  v <- pmin(pmax(v, 0), 1)
  if (isTRUE(normalize) && sum(v) > 0) v <- v / sum(v)
  
  df <- tibble(base = factor(req, levels = req), value = v)
  
  pal <- c(BASE4_NAR, "No add" = "#7F7F7F")
  
  ggplot(df, aes(base, value, fill = base)) +
    geom_col(width = 0.70, colour = "black", linewidth = 0.4) +
    scale_fill_manual(values = pal, drop = FALSE) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.12))
    ) +
    labs(x = NULL, y = "Addition probability", fill = NULL) +
    theme_nar(base_size = base_size, base_family = base_family,
              legend_inside = legend_inside, legend_pos = legend_pos) +
    theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      axis.ticks.x = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
}

library(ggrepel)

plot_pca_nar <- function(
    par_lst,
    col_lab = NULL,
    base_family = "Arial",
    base_size = 7
) {
  X3 <- sapply(par_lst, function(x) x$eff_f3)
  X5 <- X3[rownames(X3) != "AAA", , drop = FALSE]
  Xs <- scale(t(X5), center = TRUE, scale = TRUE)
  
  pc <- stats::prcomp(Xs, center = FALSE, scale. = FALSE)
  var_expl <- 100 * (pc$sdev^2 / sum(pc$sdev^2))
  
  df <- as.data.frame(pc$x[, 1:3, drop = FALSE])
  df$sample <- rownames(df)
  
  # Automatic lab (according to the rules you have always been in front of you)
  df$lab <- dplyr::case_when(
    stringr::str_detect(df$sample, stringr::regex("RNase\\s*I", ignore_case = TRUE)) ~ "RNase I",
    stringr::str_detect(df$sample, stringr::regex("\\bMNase\\b", ignore_case = TRUE)) ~ "MNase",
    stringr::str_detect(df$sample, stringr::regex("\\bP1\\b", ignore_case = TRUE))    ~ "P1",
    TRUE ~ "Other"
  )
  
  df$lab <- factor(df$lab, levels = c("RNase I", "MNase", "P1", "Other"))
  
  # Default color (you can also pass Col_lab)
  if (is.null(col_lab)) {
    col_lab <- c("RNase I" = "#0072B2", "MNase" = "#009E73", "P1" = "#D55E00", "Other" = "grey50")
  }
  
  ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = lab)) +
    ggplot2::geom_point(size = 2, stroke = 0.2, show.legend = FALSE) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = sample),
      size = 2.2,
      max.overlaps = 50,
      box.padding = 0.25,
      point.padding = 0.30,
      min.segment.length = 0,
      segment.color = "grey60",
      segment.size = 0.25,
      seed = 134,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = col_lab, name = NULL, drop = FALSE) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", var_expl[2])
    ) +
    theme_nar(base_size = base_size, base_family = base_family, legend_inside = FALSE) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}

library(patchwork)

build_fig3_plot <- function(par_lst2, par_lst3) {
  phd <- plot_hd_nar(par_lst3)
  phd <- phd + ggplot2::theme(
    strip.text.y.left = ggplot2::element_text(size = 8 - 1, hjust = 1)
  )

  par_0 <- par_lst3$`Darnell, MNase`
  pbase <- plot_mnase_base_specificity_nar(par_0, normalize = "none")

  par_5add <- par_lst2$`Darnell MNase`$prob_add5
  names(par_5add)[5] <- "No add"
  p5add <- plot_add5_prob_nar(par_5add)

  plig <- plot_pca_nar(par_lst3)

  right_row1 <- (pbase + p5add) + plot_layout(ncol = 2, widths = c(1, 1))
  right_row2 <- (plig + plot_spacer()) + plot_layout(ncol = 2, widths = c(1, 0))
  right_block <- (right_row1 / right_row2) + plot_layout(heights = c(1, 2))

  (phd | right_block) +
    plot_layout(widths = c(1, 1.0)) +
    plot_annotation(tag_levels = "A")
}

run_fig3 <- function(
    project_dir = NULL,
    input_rdata = NULL,
    output_pdf = NULL) {
  if (is.null(project_dir)) {
    if (interactive()) {
      project_dir <- readline("Please input ORF_calling directory: ")
    } else {
      stop("project_dir is required. Example: run_fig3(project_dir='.')")
    }
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  if (is.null(input_rdata)) {
    input_rdata <- file.path(project_dir, "figures", "data", "fig3_input.RData")
  }
  if (is.null(output_pdf)) {
    output_pdf <- file.path(project_dir, "figures", "Fig_3_par_summary.pdf")
  }

  env <- new.env(parent = emptyenv())
  load(input_rdata, envir = env)
  if (!exists("par_lst2", envir = env) || !exists("par_lst3", envir = env)) {
    stop("Input RData must contain par_lst2 and par_lst3.")
  }

  p <- build_fig3_plot(env$par_lst2, env$par_lst3)
  ggsave(
    filename = output_pdf,
    plot = p,
    width = 178, height = 120, units = "mm",
    device = grDevices::cairo_pdf,
    dpi = 600, bg = "white"
  )
  message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
  invisible(p)
}

if (sys.nframe() == 0) {
  run_fig3()
}
