#!/usr/bin/env Rscript

infer_project_dir <- function(project_dir = NULL) {
  if (!is.null(project_dir) && nzchar(project_dir)) {
    return(normalizePath(project_dir, mustWork = TRUE))
  }

  if (sys.nframe() > 0) {
    for (i in rev(seq_len(sys.nframe()))) {
      fr <- sys.frame(i)
      if (exists("ofile", envir = fr, inherits = FALSE)) {
        ofile <- get("ofile", envir = fr, inherits = FALSE)
        if (is.character(ofile) && length(ofile) > 0 && nzchar(ofile[1]) && ofile[1] != "-" && file.exists(ofile[1])) {
          script_path <- normalizePath(ofile[1], mustWork = TRUE)
          script_dir <- dirname(script_path)
          if (basename(script_dir) == "figures") {
            return(dirname(script_dir))
          }
        }
      }
    }
  }

  normalizePath(getwd(), mustWork = TRUE)
}

resolve_external_data_dir <- function(project_dir) {
  # Highest priority: explicit override.
  env_dir <- Sys.getenv("RIBOBA_DATA_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = TRUE))
  }

  # Otherwise probe common layouts. Each candidate is the directory that CONTAINS
  # figure_ready_data/, so the Zenodo archives work after a plain unzip with no
  # environment variable. (default_rds already carries the "figure_ready_data/"
  # prefix, so we return the parent, not figure_ready_data itself.)
  candidates <- c(
    project_dir,                                              # data placed inside code repo
    file.path(project_dir, ".."),                             # code + data unzipped side by side
    file.path(project_dir, "..", "RiboBA_analysis_data_supp"),# sibling (this supplement)
    file.path(project_dir, "..", "RiboBA_analysis_data"),     # sibling (original data repo)
    file.path(project_dir, "..", "..", "RiboBA_analysis_data_supp"),
    file.path(project_dir, "..", "..", "RiboBA_analysis_data")
  )
  for (cand in candidates) {
    p <- normalizePath(cand, mustWork = FALSE)
    if (dir.exists(file.path(p, "figure_ready_data"))) {
      return(p)
    }
  }

  ""
}

resolve_rds_path <- function(default_rds, project_dir) {
  if (grepl("^/", default_rds)) {
    return(default_rds)
  }

  local_path <- file.path(project_dir, default_rds)
  if (file.exists(local_path)) {
    return(local_path)
  }

  data_dir <- resolve_external_data_dir(project_dir)
  if (!nzchar(data_dir)) {
    return(local_path)
  }

  file.path(data_dir, default_rds)
}

find_data_list <- function(required, project_dir, default_rds = NULL, env = .GlobalEnv) {
  candidates <- c("fig_data", "plot_data", "input_data", ".Last.value")

  for (nm in candidates) {
    if (exists(nm, envir = env, inherits = TRUE)) {
      obj <- get(nm, envir = env, inherits = TRUE)
      if (is.list(obj) && all(required %in% names(obj))) return(obj)
    }
  }

  # Fallback: search all global variables for a matching named list
  for (nm in ls(envir = env, all.names = TRUE)) {
    obj <- get(nm, envir = env, inherits = FALSE)
    if (is.list(obj) && all(required %in% names(obj))) return(obj)
  }

  if (!is.null(default_rds)) {
    rds_path <- resolve_rds_path(default_rds = default_rds, project_dir = project_dir)
    if (!file.exists(rds_path)) stop("Default RDS not found: ", rds_path)
    obj <- readRDS(rds_path)
    if (!is.list(obj)) stop("RDS must contain a named list: ", rds_path)
    miss <- setdiff(required, names(obj))
    if (length(miss) > 0L) stop("RDS missing objects: ", paste(miss, collapse = ", "))
    return(obj)
  }

  stop(
    "Cannot find input data list in current session. ",
    "Please either readRDS(...) first or provide a valid default RDS mapping."
  )
}

write_list_rdata <- function(x, out_file) {
  e <- new.env(parent = emptyenv())
  list2env(x, envir = e)
  save(list = names(x), file = out_file, envir = e)
}

run_from_rds_source <- function(script_rel,
                                run_fun,
                                required = character(0),
                                default_rds = NULL,
                                output_pdf = NULL,
                                project_dir = NULL) {
  project_dir <- infer_project_dir(project_dir)
  script_path <- file.path(project_dir, script_rel)
  if (!file.exists(script_path)) stop("Script not found: ", script_path)

  source(script_path, local = .GlobalEnv)
  if (!exists(run_fun, envir = .GlobalEnv, inherits = FALSE)) {
    stop("Run function not found after sourcing script: ", run_fun)
  }
  run_obj <- get(run_fun, envir = .GlobalEnv, inherits = FALSE)

  if (length(required) == 0L) {
    args <- list(project_dir = project_dir)
    if (!is.null(output_pdf)) args$output_pdf <- file.path(project_dir, output_pdf)
    do.call(run_obj, args)
    return(invisible(TRUE))
  }

  x <- find_data_list(required = required, project_dir = project_dir, default_rds = default_rds)
  tmp_rdata <- tempfile(fileext = ".RData")
  on.exit(unlink(tmp_rdata), add = TRUE)
  write_list_rdata(x, tmp_rdata)

  args <- list(project_dir = project_dir, input_rdata = tmp_rdata)
  if (!is.null(output_pdf)) args$output_pdf <- file.path(project_dir, output_pdf)
  do.call(run_obj, args)

  invisible(TRUE)
}
