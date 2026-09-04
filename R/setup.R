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
  configure_compute(cfg, quiet = quiet)
  if (!quiet) log_step("project root: ", PROJECT_ROOT)
  cfg
}

#' Apply the compute settings from the config.
#'
#' Seurat dispatches SCTransform and friends through `future`, which caps how
#' much captured data it will hand to an evaluation at 500 MiB by default.
#' SCTransform's conserve.memory path captures the count matrix in its closure
#' and blows past that on a full sample, so the cap is raised here, once, for
#' every script. The cap is a transfer ceiling rather than an allocation:
#' raising it reserves nothing.
configure_compute <- function(cfg, quiet = FALSE) {
  compute <- cfg$compute %||% list()
  max_size_gb <- compute$globals_max_size_gb %||% 8

  options(future.globals.maxSize = max_size_gb * 1024^3)
  if (!quiet) log_step(sprintf("future.globals.maxSize = %g GiB", max_size_gb))

  plan_name <- compute$future_plan %||% "sequential"
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(plan_name)
  } else if (!identical(plan_name, "sequential")) {
    warning("future is not installed; compute.future_plan is ignored")
  }
  invisible(cfg)
}

#' Source every analysis module. Kept separate from init_project() so that
#' individual modules can be sourced on their own during development.
load_modules <- function(modules = c("io", "preprocess", "dimred", "annotate",
                                     "composition", "trajectory", "de", "gsea",
                                     "plots")) {
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

# ---------------------------------------------------------------------------
# Object helpers
#
# These exist because of a specific failure mode. A missing metadata column
# does not error in R: `obj$absent == "x"` returns logical(0), and
# `TRUE & logical(0)` is logical(0), so one absent column silently empties an
# entire multi-condition filter. subset() then reports "No cells found", which
# points at the data when the actual problem is a column an earlier step never
# wrote. select_cells() turns that into a message naming the condition.
# ---------------------------------------------------------------------------

#' Stop unless the object carries these metadata columns.
require_metadata <- function(obj, columns, context = "this step") {
  absent <- setdiff(columns, colnames(obj@meta.data))
  if (length(absent))
    stop(context, " needs metadata column(s) the object does not have: ",
         paste(absent, collapse = ", "),
         "\nThe object has: ", paste(colnames(obj@meta.data), collapse = ", "),
         "\nAn earlier pipeline step should have written them; re-run it and ",
         "check its output for warnings.")
  invisible(TRUE)
}

#' Subset to cells passing every named condition, reporting the breakdown.
#'
#' @param conditions named list of logical vectors, one value per cell. NA
#'   counts as a failure. The names are what gets reported, so make them read
#'   like the condition they describe.
select_cells <- function(obj, conditions, context = "subset", quiet = FALSE) {
  n <- ncol(obj)

  wrong_length <- vapply(conditions, function(x) length(x) != n, logical(1))
  if (any(wrong_length))
    stop(context, ": condition(s) [",
         paste(names(conditions)[wrong_length], collapse = ", "),
         "] produced ",
         paste(vapply(conditions[wrong_length], length, integer(1)), collapse = "/"),
         " values for ", n, " cells.",
         "\nA metadata column that does not exist yields a zero-length ",
         "condition, which would otherwise empty the filter silently. ",
         "Check that the column names in this filter exist on the object.")

  passes <- lapply(conditions, function(x) !is.na(x) & x)

  # Cumulative counts show which condition does the damage, which a per-
  # condition count alone does not.
  cumulative <- Reduce(`&`, passes, accumulate = TRUE)
  report <- data.frame(
    condition = names(conditions),
    passing = vapply(passes, sum, integer(1)),
    remaining = vapply(cumulative, sum, integer(1)),
    row.names = NULL
  )

  keep <- cumulative[[length(cumulative)]]
  if (!any(keep)) {
    stop(context, ": no cells passed all conditions (starting from ", n, ").\n",
         paste(utils::capture.output(print(report, row.names = FALSE)),
               collapse = "\n"),
         "\nThe first row where `remaining` hits 0 is the condition to look at.")
  }

  if (!quiet) {
    log_step(context, ": ", sum(keep), " of ", n, " cells")
    print(report, row.names = FALSE)
  }
  subset(obj, cells = colnames(obj)[keep])
}

#' How many cells carry a non-NA value for each column.
metadata_coverage <- function(obj, columns) {
  present <- intersect(columns, colnames(obj@meta.data))
  data.frame(
    column = columns,
    present = columns %in% present,
    non_na = vapply(columns, function(col) {
      if (!col %in% present) return(NA_integer_)
      sum(!is.na(obj@meta.data[[col]]))
    }, integer(1)),
    row.names = NULL
  )
}
