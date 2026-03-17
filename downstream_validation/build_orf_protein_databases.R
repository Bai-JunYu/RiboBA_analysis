#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

source(file.path("figures", "helpers", "prepare_utils.R"))

load_tool_lists <- function(public_inputs_dir) {
  e <- new.env(parent = emptyenv())
  load(file.path(public_inputs_dir, "result", "p3", "for_orfrater", "ribo_orfrater.RData"), envir = e)
  load(file.path(public_inputs_dir, "result", "p3", "price", "ribo_price.RData"), envir = e)
  load(file.path(public_inputs_dir, "result", "p3", "ribocode", "ribo_ribocode.RData"), envir = e)
  load(file.path(public_inputs_dir, "result", "p3", "ribotish", "ribo_ribotish.RData"), envir = e)
  load(file.path(public_inputs_dir, "result", "p3", "ribORF", "ribo_riborf.RData"), envir = e)
  load(file.path(public_inputs_dir, "data", "p1", "riboba_orf.RData"), envir = e)

  orf_riboba <- lapply(e$riboba_orf, function(x) x$orf_info$ORF)

  list(
    PRICE = e$orf_price,
    RiboBA = orf_riboba,
    RiboCode = e$orf_ribocode,
    RiboTISH = e$orf_ribotish,
    RibORF = e$orf_riborf,
    `ORF-RATER` = e$orf_rater
  )
}

build_df_ms_cnt <- function(tool_lst, public_inputs_dir) {
  root <- file.path(public_inputs_dir, "data", "p3", "HLA_ms")
  srr_set <- c("SRR23242344", "SRR23242345", "SRR23242346", "SRR23242347", "SRR2433794", "SRR7073124", "SRR8449566", "SRR8449567", "SRR8449568")
  tool_dirs <- c("riboba", "price", "ribocode", "ribotish", "riborf", "rater")
  tool_lv <- c("RiboBA", "RiboTISH", "ORF-RATER", "RiboCode", "PRICE", "RibORF")

  canon_tool <- function(x) dplyr::recode(tolower(x), riboba = "RiboBA", price = "PRICE", ribocode = "RiboCode", ribotish = "RiboTISH", riborf = "RibORF", rater = "ORF-RATER", .default = x)

  parse_tool_srr <- function(p) {
    comps <- strsplit(normalizePath(p, mustWork = FALSE), "/")[[1]]
    toolseg <- comps[tolower(comps) %in% tool_dirs]
    toolseg <- if (length(toolseg)) toolseg[1] else NA_character_
    srr_all <- regmatches(p, gregexpr("SRR\\d+", p))[[1]]
    srrseg <- if (length(srr_all)) tail(srr_all, 1) else NA_character_
    list(tool = canon_tool(toolseg), srr = srrseg)
  }

  read_one_msstats <- function(p) {
    meta <- parse_tool_srr(p)
    tibble(tool = meta$tool, srr = meta$srr) %>% bind_cols(read.csv(p, sep = ","))
  }

  ms_paths <- list.files(root, pattern = "msstats\\.csv$", recursive = TRUE, full.names = TRUE)
  ms_paths <- ms_paths[grepl(paste(tool_dirs, collapse = "|"), ms_paths, ignore.case = TRUE) & grepl(paste(srr_set, collapse = "|"), ms_paths, ignore.case = TRUE)]
  if (!length(ms_paths)) stop("No msstats.csv found under: ", root)

  types_keep <- c("uORF", "uoORF", "ncRNA", "intORF", "doORF", "dORF")

  ms_valid <- bind_rows(lapply(ms_paths, read_one_msstats)) %>%
    filter(Protein.Description %in% c("true", TRUE)) %>%
    filter(grepl(paste(types_keep, collapse = "|"), Protein)) %>%
    mutate(orf_id_str = sub("^[^_]+_([^_]+)_.*$", "\\1", Protein)) %>%
    filter(srr %in% srr_set) %>%
    distinct(tool, srr, orf_id_str)

  keep_cols <- c("orf_biotype", "orf_id", "is_inframe")
  id_col_guess <- c("orf_id", "orf_id_str", "ORF_id", "orfID", "id")

  pull_rows_by_id <- function(df, ids, id_cols = id_col_guess) {
    if (is.null(df) || !nrow(df)) return(df[0, , drop = FALSE])
    hit <- intersect(id_cols, names(df))
    if (!length(hit)) return(df[0, , drop = FALSE])
    id_col <- hit[1]
    df[as.character(df[[id_col]]) %in% as.character(ids), , drop = FALSE]
  }

  supp_ms_orf <- ms_valid %>%
    mutate(tool = as.character(tool), srr = as.character(srr), orf_id_str = as.character(orf_id_str)) %>%
    filter(tool %in% names(tool_lst)) %>%
    group_by(tool, srr) %>%
    group_modify(~{
      df_tool <- tool_lst[[.y$tool]][[.y$srr]]
      ids <- .x$orf_id_str
      if (is.null(df_tool) || !nrow(df_tool)) return(tibble(orf_id_str = ids, matched = FALSE))
      got <- pull_rows_by_id(df_tool, ids)
      if (!nrow(got)) return(tibble(orf_id_str = ids, matched = FALSE))
      got <- as_tibble(got)
      hit <- intersect(id_col_guess, names(got))
      df_idcol <- if (length(hit)) hit[1] else NA_character_
      if (!is.na(df_idcol)) got <- got %>% mutate(orf_id_str = as.character(.data[[df_idcol]]))
      keep_here <- intersect(keep_cols, names(got))
      got %>% select(any_of(keep_here), orf_id_str) %>% mutate(matched = TRUE)
    }) %>%
    ungroup() %>%
    mutate(tool = factor(tool, levels = tool_lv))

  biotype_norf <- c("uORF", "uoORF", "ncRNA", "intORF", "dORF", "doORF")
  biotype_other <- c("CDS", "Ext", "Trunc", "Variant")

  supp_ms_orf_keep <- supp_ms_orf %>%
    filter((orf_biotype %in% biotype_norf & is_inframe == FALSE) | (orf_biotype %in% biotype_other))

  biotype_lv <- c("uORF", "uoORF", "ncRNA", "intORF", "dORF", "doORF")
  supp_ms_orf_keep %>%
    filter(orf_biotype %in% biotype_lv) %>%
    count(tool, orf_biotype, name = "n") %>%
    mutate(tool = factor(tool, levels = tool_lv), orf_biotype = factor(orf_biotype, levels = biotype_lv))
}

run_prepare_fig_s8_data <- function(project_dir = NULL, workspace_rdata = NULL) {
  project_dir <- get_project_dir(project_dir)
  public_inputs_dir <- resolve_public_inputs_dir(project_dir)

  tool_lst <- load_tool_lists(public_inputs_dir)
  df_ms_cnt <- build_df_ms_cnt(tool_lst, public_inputs_dir)

  out <- file.path(project_dir, "figures", "data", "fig_s8_input.RData")
  save(df_ms_cnt, file = out)
  message("Saved: ", normalizePath(out, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  args <- parse_prepare_args()
  run_prepare_fig_s8_data(args$project_dir, args$workspace_rdata)
}
