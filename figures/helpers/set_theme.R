# ===== style kit =====
library(ggplot2)
library(grid)
library(svglite)

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

  # `source()` stores the source file path in a frame-local `ofile`.
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

SCRIPT_DIR <- get_script_dir()
PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
PUBLIC_INPUTS_DIR <- Sys.getenv("PUBLIC_INPUTS_DIR", unset = file.path(PROJECT_DIR, "public_inputs"))
RIBOBC_DIR <- Sys.getenv("RIBOBC_DIR", unset = file.path(PROJECT_DIR, "RiboBC"))
RIBOCONTEXT_DIR <- Sys.getenv("RIBOCONTEXT_DIR", unset = file.path(PROJECT_DIR, "RiboContext"))

# 1) Unified Theme
theme_nar <- function(
    base_size = 8,
    base_family = "Arial",
    line_w = 0.3,
    grid_w = 0.25,
    grid_col = "grey88",
    legend_inside = TRUE,
    legend_pos = c(0.02, 0.98),
    legend_bg_alpha = 0.85
) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      # Global force font (avoid patchwork/annotation fallback to sans)
      text = ggplot2::element_text(family = base_family),
      
      # Grid: Keep major only
      panel.grid.major = ggplot2::element_line(linewidth = grid_w, colour = grid_col),
      panel.grid.minor = ggplot2::element_blank(),
      
      # Panels & Sections
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain", size = base_size),
      panel.spacing = grid::unit(3, "mm"),
      
      # Key: Avoid left/bottom overlay thickening - only draw panel.border, close axis.line
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = line_w),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(linewidth = line_w, colour = "black"),
      axis.ticks.length = grid::unit(1.1, "mm"),
      
      axis.title = ggplot2::element_text(size = base_size),
      axis.text  = ggplot2::element_text(size = base_size - 1),
      
      # Legend: default inside (you can also legend_inside = false)
      legend.position = if (legend_inside) "inside" else "right",
      legend.position.inside = legend_pos,
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(size = base_size),
      legend.text  = ggplot2::element_text(size = base_size - 1, lineheight = 1.12),
      legend.key.width  = grid::unit(10, "pt"),
      legend.key.height = grid::unit(9, "pt"),
      legend.spacing.y  = grid::unit(1.5, "pt"),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", legend_bg_alpha),
        colour = NA, linewidth = 0
      ),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA),
      
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
}

# 2) Uniform color matching (color blindness friendly, default Okabe-Ito, can be changed by yourself)
nar_palette <- function() {
  c(
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#009E73", # bluish green
    "#CC79A7", # reddish purple
    "#F0E442", # yellow
    "#56B4E9", # sky blue
    "#000000", # black
    "#E69F00"  # orange
  )
}

scale_color_nar <- function(values = NULL, ...) {
  if (is.null(values)) values <- nar_palette()
  ggplot2::scale_colour_manual(values = values, ...)
}

scale_fill_nar <- function(values = NULL, ...) {
  if (is.null(values)) values <- nar_palette()
  ggplot2::scale_fill_manual(values = values, ...)
}
# 3) Save uniformly (size mm; PDF force embedded font)
save_nar <- function(filename, plot,
                     width_mm = 178, height_mm = 230,
                     dpi = 600,
                     family = "Arial") {
  ggplot2::ggsave(
    filename, plot,
    width = width_mm, height = height_mm, units = "mm",
    device = svglite,
    family = family,
    dpi = dpi
  )
}
