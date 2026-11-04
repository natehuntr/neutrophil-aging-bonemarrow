#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Build analysis/pipeline.qmd: the whole pipeline as ONE Quarto file.
#
#   Rscript scripts/build_pipeline_qmd.R
#
# The single file is generated, not maintained. Every R/ module and every
# step script is copied in verbatim, with only the lines that locate files on
# disk replaced (the project root, the config file, `source()` calls), so the
# one-file version runs exactly the code the multi-file pipeline runs. Re-run
# this after changing anything in R/, scripts/ or config/ and commit the
# result alongside.
#
# Writes: analysis/pipeline.qmd
# ---------------------------------------------------------------------------

root <- if (file.exists("config/config.yml")) "." else ".."
read <- function(...) readLines(file.path(root, ...), warn = FALSE)

#' Replace an exact block of text, failing loudly if the source has moved on.
replace_once <- function(lines, old, new, what) {
  text <- paste(lines, collapse = "\n")
  hits <- gregexpr(old, text, fixed = TRUE)[[1]]
  if (length(hits) != 1 || hits[1] == -1)
    stop("build_pipeline_qmd: expected exactly one '", what, "' in the source; ",
         "the code changed shape -- update this generator")
  strsplit(sub(old, new, text, fixed = TRUE), "\n", fixed = TRUE)[[1]]
}

#' The lines every script uses to set itself up as a standalone file.
strip_standalone_setup <- function(lines) {
  drop <- grepl('^source\\(if \\(file.exists\\("R/setup.R"\\)\\)', lines) |
    grepl("^cfg <- init_project\\(.*\\)\\s*$", lines) |
    grepl("^load_modules\\(\\)\\s*$", lines) |
    grepl("^#!/usr/bin/env Rscript", lines) |
    grepl("^quit\\(status", lines)
  # Step 13 sources its plotting modules itself; they are already defined.
  src <- grep('source\\(project_path\\("R"', lines)
  for (i in src) {
    drop[i] <- TRUE
    if (i > 1 && grepl("^for \\(m in c\\(", lines[i - 1])) drop[i - 1] <- TRUE
  }
  lines[!drop]
}

chunk <- function(label, code, options = character()) {
  c(paste0("```{r ", label, "}"), paste0("#| ", options), code, "```", "")
}

# --- setup.R, adapted to a single file ---------------------------------------
setup <- read("R", "setup.R")
start <- grep("^PROJECT_ROOT <- local\\(\\{", setup)
end <- start - 1 + which(setup[start:length(setup)] == "})")[1]
setup <- c(setup[seq_len(start - 1)],
           "# Single-file mode: knitr runs every chunk in the directory holding",
           "# pipeline.qmd, so that directory is the project root.",
           "PROJECT_ROOT <- normalizePath(getwd())",
           setup[(end + 1):length(setup)])
setup <- replace_once(setup, "  cfg <- yaml::read_yaml(path)",
                      paste0("  # Single-file mode: the config is embedded above, and the render\n",
                             "  # parameters override it.\n",
                             "  cfg <- apply_render_params(yaml::yaml.load(PIPELINE_CONFIG_YAML))"),
                      "read_yaml in load_config")
setup <- replace_once(setup, '  source(project_path("R", "packages.R"))\n', "",
                      "source of packages.R in init_project")

modules <- setdiff(sub("\\.R$", "", list.files(file.path(root, "R"), "\\.R$")),
                   c("setup", "packages"))

steps <- list(
  list(n = 1,  file = "01_preprocess.R",          title = "Preprocess: ADT/HTO, demultiplexing, QC, doublets", heavy = TRUE),
  list(n = 2,  file = "02_annotate.R",            title = "Reductions, SingleR labels, CytoTRACE2", heavy = TRUE),
  list(n = 3,  file = "03_merge.R",               title = "Merge sexes, stage granulocytes, depth-matched assays", heavy = TRUE),
  list(n = 4,  file = "04_sex_candidates.R",      title = "Sex-difference candidates for orthogonal validation", heavy = TRUE),
  list(n = 5,  file = "05_trajectory.R",          title = "Trajectories and tradeSeq GAMs", heavy = TRUE),
  list(n = 6,  file = "06_age_trends.R",          title = "Age trends and pairwise age changes, RNA and ADT", heavy = TRUE),
  list(n = 7,  file = "07_glm_models.R",          title = "NB-GLM omnibus, interaction and trajectory models", heavy = TRUE),
  list(n = 8,  file = "08_gsea.R",                title = "Sex x age interaction GSEA, per stage", heavy = TRUE),
  list(n = 9,  file = "09_development_shifts.R",  title = "Stage composition, pseudotime and potency shifts", heavy = TRUE),
  list(n = 10, file = "10_endpoint_contrast.R",   title = "3m vs 18m genes and gene sets, per stage", heavy = TRUE),
  list(n = 11, file = "11_progenitors.R",         title = "The progenitor compartment", heavy = TRUE),
  list(n = 12, file = "12_protein_gate.R",        title = "Stage composition from surface protein, before and after QC", heavy = TRUE),
  list(n = 13, file = "13_figures.R",             title = "Summary figures from the tables", heavy = FALSE),
  list(n = 14, file = "collect_results.R",        title = "Results bundle", heavy = FALSE)
)

config_yaml <- read("config", "config.yml")
if (any(grepl(")-----", config_yaml, fixed = TRUE)))
  stop("config.yml contains the raw-string delimiter; pick another one here")

# --- the document ------------------------------------------------------------
out <- c(
'---',
'title: "Neutrophil ageing in mouse bone marrow"',
'subtitle: "The whole pipeline in one file: Cell Ranger output to figures"',
'format:',
'  html:',
'    toc: true',
'    toc-depth: 2',
'    code-fold: true',
'    code-summary: "Code"',
'    embed-resources: true',
'execute:',
'  warning: false',
'  message: false',
'params:',
'  # Folder holding one Cell Ranger per-sample output directory per sex',
'  # (relative to this file, or absolute).',
'  cellranger_root: "data/per_sample_outs"',
'  # Which steps to run: "1-13", "4,6,9", "13" ... Step 14 writes the bundle.',
'  steps: "1-14"',
'  # Re-run steps that already completed (each writes results/logs/stepNN.done).',
'  force: false',
'  # Cores for the packages that fork workers (CytoTRACE2, monocle3).',
'  cores: 4',
'  # Install the R packages first (CRAN, Bioconductor and GitHub).',
'  install: false',
'---',
'',
'<!-- GENERATED by scripts/build_pipeline_qmd.R from R/, scripts/ and',
'     config/config.yml. Edit those and re-run the generator; edits made here',
'     are lost the next time it runs. -->',
'',
'## How to run this',
'',
'This file is the complete analysis: every function and every step of the',
'pipeline, reading the Cell Ranger output and writing the same tables, objects',
'and figures as the multi-file version in the repository. Put it in an empty',
'folder (or leave it in `analysis/`), point it at the data, and render:',
'',
'```bash',
'quarto render pipeline.qmd -P cellranger_root:/path/to/per_sample_outs',
'quarto render pipeline.qmd -P steps:4-14            # later steps only',
'quarto render pipeline.qmd -P steps:13              # redraw the figures',
'quarto render pipeline.qmd -P install:true -P steps:none   # install packages only',
'```',
'',
'Everything is written next to this file, under `results/` (`objects/`,',
'`tables/`, `figures/`, `logs/`, `summary/`).',
'',
'**Steps are resumable.** A step that completes writes',
'`results/logs/stepNN.done` and is skipped on the next render, so a render that',
'stops part-way can simply be started again. Set `force: true` to re-run',
'completed steps. A step that fails is reported in this document with the',
'end of its log, and the render carries on; steps that need its output then',
'fail in turn, and say which object they are missing.',
'',
'**Resources.** Steps 1-3 and 5 need the most memory (plan on 64 GB for',
'the full merged object), and step 7 is the slowest by far: it refits every',
'model `glm_de.n_perm` times (20 in the config below). On a laptop, lower',
'`glm_de.n_perm`, `null_model.n_splits` and `endpoint.n_perm` in the',
'configuration chunk for a first pass. Step 13 needs no objects and runs in',
'seconds.',
'',
'**Requirements.** R >= 4.3 with Seurat v5, and the packages listed in the',
'installation chunk (`install: true` installs them). Leiden clustering',
'(`clustering.algorithm: 4`) also needs the Python `leidenalg` package through',
'reticulate; set it to 1 for Louvain.',
'',
'## Configuration',
'',
'The project configuration, embedded. Every path, threshold and tuning',
'parameter of the analysis is here; the render parameters above override',
'`paths.cellranger_root` and `compute.cores`.',
'',
chunk("config", c(
  "PIPELINE_CONFIG_YAML <- r\"-----(",
  config_yaml,
  ")-----\"",
  "",
  "#' Apply the render parameters on top of the embedded configuration.",
  "apply_render_params <- function(cfg) {",
  "  p <- if (exists(\"params\")) params else list()",
  "  if (!is.null(p$cellranger_root) && nzchar(p$cellranger_root))",
  "    cfg$paths$cellranger_root <- p$cellranger_root",
  "  if (!is.null(p$cores)) cfg$compute$cores <- as.integer(p$cores)",
  "  cfg",
  "}")),
'## Package installation',
'',
'Runs only with `install: true`. Installs CRAN, Bioconductor and GitHub',
'packages into the first writable library (`R_LIBS_USER` if set). Set',
'`GITHUB_PAT` in the environment to avoid GitHub rate limits.',
'',
chunk("install", strip_standalone_setup(read("scripts", "00_install_dependencies.R")),
      c("eval: !expr isTRUE(params$install)", "code-fold: true")),
'## Functions',
'',
'Project setup, then one chunk per module of the repository\'s `R/` folder.',
'These only define functions and constants; nothing runs until the steps.',
'',
chunk("module-setup", setup),
chunk("module-packages", gsub("Run: Rscript scripts/00_install_dependencies.R",
                               "Render with -P install:true to install them.",
                               read("R", "packages.R"), fixed = TRUE)),
unlist(lapply(modules, function(m) chunk(paste0("module-", m), read("R", paste0(m, ".R"))))),
'## Running the steps',
'',
chunk("step-runner", c(
  "cfg <- load_config()",
  "ensure_output_dirs(cfg)",
  "LOG_DIR <- project_path(\"results\", \"logs\")",
  "dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)",
  "",
  "#' Which steps the `steps` parameter asks for: \"1-13\", \"4,6,9\", \"none\".",
  "parse_steps <- function(spec) {",
  "  spec <- gsub(\"\\\\s\", \"\", as.character(spec))",
  "  if (!nzchar(spec) || identical(spec, \"none\")) return(integer())",
  "  unlist(lapply(strsplit(spec, \",\")[[1]], function(part) {",
  "    ends <- as.integer(strsplit(part, \"-\")[[1]])",
  "    if (length(ends) == 2) seq(ends[1], ends[2]) else ends",
  "  }))",
  "}",
  "WANTED <- parse_steps(params$steps)",
  "STATUS <- list()",
  "",
  "#' Run one step in its own environment, logging to results/logs/stepNN.log.",
  "#'",
  "#' Each step's code is the repository script, unchanged apart from its",
  "#' setup lines. It runs in a fresh environment so one step's variables",
  "#' cannot leak into the next, as when each script ran in its own process.",
  "run_step <- function(n, title, code, needs_packages = TRUE) {",
  "  id <- sprintf(\"step%02d\", n)",
  "  done <- file.path(LOG_DIR, paste0(id, \".done\"))",
  "  log <- file.path(LOG_DIR, paste0(id, \".log\"))",
  "  record <- function(status, mins = NA_real_)",
  "    STATUS[[id]] <<- data.frame(step = n, title = title, status = status,",
  "                                minutes = round(mins, 1), stringsAsFactors = FALSE)",
  "  if (!n %in% WANTED) { record(\"not requested\"); return(invisible()) }",
  "  if (file.exists(done) && !isTRUE(params$force)) {",
  "    cat(\"Already completed (\", format(file.mtime(done)), \"); skipped. \",",
  "        \"Render with force: true to re-run.\\n\", sep = \"\")",
  "    record(\"skipped: already done\"); return(invisible())",
  "  }",
  "  set.seed(cfg$seed)",
  "",
  "  env <- new.env(parent = globalenv())",
  "  # Scripts read their command-line arguments; a render has none to give.",
  "  env$commandArgs <- function(...) character(0)",
  "  environment(code) <- env",
  "",
  "  con <- file(log, open = \"wt\")",
  "  # Printed output goes through the sink. Messages and warnings are caught",
  "  # by a handler instead: during a render knitr intercepts them before any",
  "  # message sink would see them, and the log came out empty.",
  "  sink(con)",
  "  to_log <- function(prefix) function(cond) {",
  "    cat(prefix, conditionMessage(cond), if (prefix != \"\") \"\\n\", file = con, sep = \"\")",
  "    invokeRestart(if (inherits(cond, \"warning\")) \"muffleWarning\" else \"muffleMessage\")",
  "  }",
  "  started <- Sys.time()",
  "  # Package loading is inside the handler: a missing package fails this",
  "  # step, with the install hint, rather than stopping the whole render.",
  "  err <- tryCatch(withCallingHandlers({",
  "    if (needs_packages) load_packages()",
  "    configure_compute(cfg, quiet = TRUE)",
  "    code()",
  "    NULL",
  "  }, message = to_log(\"\"), warning = to_log(\"Warning: \")),",
  "  error = function(e) e)",
  "  if (!is.null(err)) cat(\"ERROR: \", conditionMessage(err), \"\\n\", file = con, sep = \"\")",
  "  sink(); close(con)",
  "  mins <- as.numeric(difftime(Sys.time(), started, units = \"mins\"))",
  "",
  "  lines <- readLines(log, warn = FALSE)",
  "  if (is.null(err)) {",
  "    writeLines(format(Sys.time()), done)",
  "    record(\"completed\", mins)",
  "    cat(sprintf(\"Completed in %.1f min. Last lines of %s:\\n\\n\", mins,",
  "                sub(paste0(PROJECT_ROOT, \"/\"), \"\", log, fixed = TRUE)))",
  "    cat(utils::tail(lines, 25), sep = \"\\n\")",
  "  } else {",
  "    unlink(done)",
  "    record(paste(\"FAILED:\", strsplit(conditionMessage(err), \"\\n\")[[1]][1]), mins)",
  "    cat(sprintf(\"FAILED after %.1f min: %s\\n\\nLast lines of the log:\\n\\n\",",
  "                mins, conditionMessage(err)))",
  "    cat(utils::tail(lines, 40), sep = \"\\n\")",
  "  }",
  "  invisible()",
  "}")),
unlist(lapply(steps, function(s) {
  body <- strip_standalone_setup(read("scripts", s$file))
  c(sprintf("### Step %d: %s", s$n, s$title), "",
    chunk(sprintf("step-%02d", s$n),
          c(sprintf("run_step(%d, \"%s\", needs_packages = %s, function() {",
                    s$n, s$title, if (s$heavy) "TRUE" else "FALSE"),
            body, "})"),
          c("comment: \"\"")))
})),
'### Status',
'',
chunk("status", c(
  "status <- do.call(rbind, STATUS)",
  "knitr::kable(status, row.names = FALSE)"),
  "echo: false"),
'## Results',
'',
'The summary figures step 13 drew, in the order a reader needs them. Each is',
'drawn from the results table named in `results/figures/summary/FIGURES.md`.',
'',
'**What the design supports.** Each sex is one 10x library, so every',
'male-versus-female number is also a library-versus-library number. Ages are',
'hashed within a library, so age comparisons within a sex are free of batch.',
'Each age is one hashtag pooling six mice: intervals describe sampling of',
'cells, not variation between mice. The results are hypotheses for per-mouse',
'validation.',
'',
chunk("results-figures", c(
  "fig_dir <- file.path(cfg$paths$figures, \"summary\")",
  "show_fig <- function(id, text) {",
  "  f <- list.files(fig_dir, paste0(\"^\", id, \"_.*\\\\.png$\"), full.names = TRUE)",
  "  cat(\"\\n\\n\", text, \"\\n\\n\", sep = \"\")",
  "  if (length(f)) cat(sprintf(\"![](%s)\\n\\n\", f[1]))",
  "  else cat(sprintf(\"> *%s not drawn: its table has not been produced.*\\n\\n\", id))",
  "}",
  "figs <- list(",
  "  fig01 = \"**Depth by stage, age and library**, on the counts the analyses read.\",",
  "  fig02 = \"**Depth matching per stratum**: each thinned toward its own shallowest group.\",",
  "  fig03 = \"**RNA against protein stage.** Disagreement concentrates at immature/mature.\",",
  "  fig04 = \"**Stage mix by age.** Ages share a library; the two panels are two libraries.\",",
  "  fig05 = \"**Stage shares across age**, with the log-odds change per age step.\",",
  "  fig06 = \"**Stages over- and under-represented by age**, within each library.\",",
  "  fig07 = \"**Progenitor compartment by age**: HSC expansion and MDP loss at 18m.\",",
  "  fig08 = \"**Progenitor shares across age.**\",",
  "  fig09 = \"**Stage mix from surface protein, before and after RNA QC.**\",",
  "  fig10 = \"**How much of each protein-defined stage RNA QC kept, by age.** A fall in one library means QC, not the marrow, removed those cells.\",",
  "  fig11 = \"**Pseudotime relative to 3m**, pooled and within stage.\",",
  "  fig12 = \"**Potency relative to 3m**, pooled and within stage.\",",
  "  fig13 = \"**Genes changing 3m to 18m per stratum**, with the detection check.\",",
  "  fig14 = \"**The 3m vs 18m genes through every age.**\",",
  "  fig15 = \"**Gene sets moving 3m to 18m per stratum.**\",",
  "  fig16 = \"**Four-age trends above the permutation null.**\",",
  "  fig17 = \"**Surface proteins changing with age.**\",",
  "  fig18 = \"**Genes changing with age, female** (stage-adjusted NB-GLM).\",",
  "  fig19 = \"**Genes changing with age, male** (stage-adjusted NB-GLM).\",",
  "  fig20 = \"**Developmental profiles changing with age or sex** (trajectory model).\",",
  "  fig21 = \"**Pathway trends differing between the libraries**, per stage.\",",
  "  fig22 = \"**Sex-difference candidates** for orthogonal validation. Not a test of sex.\",",
  "  fig23 = \"**The sex contrast against random splits of one library.**\")",
  "for (id in names(figs)) show_fig(id, figs[[id]])"),
  c("echo: false", "results: asis"))
)

dest <- file.path(root, "analysis", "pipeline.qmd")
writeLines(out, dest)
cat("wrote", dest, "(", length(out), "lines )\n")
