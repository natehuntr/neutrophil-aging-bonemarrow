# ---------------------------------------------------------------------------
# Project setup: configuration, paths, logging.
#
# Every script starts with `source("R/setup.R")`, which is the only file that
# needs to know where the project root is.
# ---------------------------------------------------------------------------

PROJECT_ROOT <- local({
  # here::here() if available, otherwise walk up for the config file so the
  # scripts work from Rscript, RStudio and Quarto alike.
  if (requireNamespace("here", quietly = TRUE)) return(here::here())
  dir <- normalizePath(getwd(), mustWork = TRUE)
  while (!file.exists(file.path(dir, "config", "config.yml"))) {
    parent <- dirname(dir)
    if (identical(parent, dir))
      stop("could not locate the project root (no config/config.yml above ", getwd(), ")")
    dir <- parent
  }
  dir
})

project_path <- function(...) file.path(PROJECT_ROOT, ...)

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Read config/config.yml and resolve its paths against the project root.
load_config <- function(path = project_path("config", "config.yml")) {
  cfg <- yaml::read_yaml(path)
  # Only the directory entries are resolved: raw_matrix_dir and
  # filtered_matrix_dir are name fragments joined onto a sample directory,
  # not paths in their own right.
  for (key in c("cellranger_root", "objects", "tables", "figures")) {
    p <- cfg$paths[[key]]
    # Absolute paths (e.g. an external Cell Ranger drive) are left alone.
    if (!grepl("^(/|[A-Za-z]:)", p)) cfg$paths[[key]] <- project_path(p)
  }
  cfg
}

#' Create the results directories named in the config.
ensure_output_dirs <- function(cfg) {
  for (p in cfg$paths[c("objects", "tables", "figures")])
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(cfg)
}

#' Standard entry point for a script: load config, seed the RNG, make dirs.
init_project <- function(quiet = FALSE) {
  source(project_path("R", "packages.R"))
  load_packages()
  cfg <- load_config()
  ensure_output_dirs(cfg)
  set.seed(cfg$seed)
  if (!quiet) log_step("project root: ", PROJECT_ROOT)
  cfg
}

#' Source every analysis module. Kept separate from init_project() so that
#' individual modules can be sourced on their own during development.
load_modules <- function(modules = c("io", "preprocess", "dimred", "annotate",
                                     "trajectory", "de", "gsea", "plots")) {
  for (m in modules) source(project_path("R", paste0(m, ".R")))
  invisible(modules)
}

log_step <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

#' Path helpers, so no script ever builds an output path by hand.
object_path <- function(cfg, name) file.path(cfg$paths$objects, name)
table_path  <- function(cfg, name) file.path(cfg$paths$tables, name)
figure_path <- function(cfg, name) file.path(cfg$paths$figures, name)

#' Save a Seurat/other object and report where it went.
save_object <- function(obj, cfg, name) {
  path <- object_path(cfg, name)
  saveRDS(obj, path)
  log_step("wrote ", path)
  invisible(path)
}

read_object <- function(cfg, name) {
  path <- object_path(cfg, name)
  if (!file.exists(path))
    stop("missing ", path, "\nRun the earlier pipeline step that creates it.")
  readRDS(path)
}

#' Write a data frame to results/tables and report where it went.
write_table <- function(df, cfg, name, row.names = FALSE) {
  path <- table_path(cfg, name)
  utils::write.csv(df, path, row.names = row.names)
  log_step("wrote ", path, " (", nrow(df), " rows)")
  invisible(path)
}

save_figure <- function(plot, cfg, name, width = 9, height = 7) {
  path <- figure_path(cfg, name)
  ggplot2::ggsave(path, plot, width = width, height = height)
  log_step("wrote ", path)
  invisible(path)
}

#' Fetch a sample's block from the config by key ("male" / "female").
sample_config <- function(cfg, key) {
  s <- cfg$samples[[key]]
  if (is.null(s)) stop("no sample '", key, "' in config.yml")
  s$key <- key
  s
}
