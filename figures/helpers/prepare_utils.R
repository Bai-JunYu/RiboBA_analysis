get_project_dir <- function(project_dir = NULL) {
  if (is.null(project_dir) || !nzchar(project_dir)) {
    if (interactive()) {
      project_dir <- readline("Please input project directory: ")
    } else {
      stop("project_dir is required")
    }
  }
  normalizePath(project_dir, mustWork = TRUE)
}

resolve_public_inputs_dir <- function(project_dir, public_inputs_dir = NULL) {
  candidates <- c(
    public_inputs_dir,
    Sys.getenv("PUBLIC_INPUTS_DIR", unset = ""),
    file.path(project_dir, "..", "public_inputs"),
    file.path(project_dir, "..", "..", "public_inputs"),
    file.path(project_dir, "public_inputs"),
    file.path(dirname(project_dir), "public_inputs")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[dir.exists(candidates)][1]
  if (is.na(hit)) {
    stop(
      "public_inputs directory not found.\n",
      "Set PUBLIC_INPUTS_DIR or place public_inputs near this project.\n",
      "Checked:\n- ",
      paste(candidates, collapse = "\n- ")
    )
  }
  normalizePath(hit, mustWork = TRUE)
}

resolve_workspace <- function(workspace_rdata = NULL, project_dir = NULL) {
  if (!is.null(workspace_rdata) && nzchar(workspace_rdata)) {
    return(normalizePath(workspace_rdata, mustWork = TRUE))
  }
  candidate <- file.path(project_dir, "figures", "data", "total_intermediates.RData")
  if (!file.exists(candidate)) {
    stop(
      "workspace_rdata is required for this figure.\n",
      "Provide workspace path explicitly, e.g. run_prepare_xxx(..., workspace_rdata='./workspace.RData').\n",
      "Default checked path not found: ", candidate
    )
  }
  normalizePath(candidate, mustWork = TRUE)
}

load_env <- function(rdata_path) {
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  e
}

save_from_env <- function(env, objs, output_rdata) {
  miss <- objs[!vapply(objs, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) {
    stop("Missing objects: ", paste(miss, collapse = ", "))
  }
  vals <- mget(objs, envir = env, inherits = FALSE)
  list2env(vals, envir = environment())
  save(list = objs, file = output_rdata)
  message("Saved: ", normalizePath(output_rdata, mustWork = FALSE))
  invisible(output_rdata)
}

parse_prepare_args <- function(default_project = getwd()) {
  args <- commandArgs(trailingOnly = TRUE)
  project_dir <- if (length(args) >= 1) args[[1]] else default_project
  workspace_rdata <- if (length(args) >= 2) args[[2]] else NULL
  list(project_dir = project_dir, workspace_rdata = workspace_rdata)
}
