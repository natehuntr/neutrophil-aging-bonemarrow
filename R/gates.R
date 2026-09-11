# ---------------------------------------------------------------------------
# Blocking control gates.
#
# Each gate corresponds to a way this dataset has already been observed to
# mislead. They are checked before the analyses that depend on them and are
# fatal by default, because the failure mode they guard against is a result
# that looks fine.
# ---------------------------------------------------------------------------

#' Record one gate's outcome.
gate_result <- function(name, passed, detail, blocks = "") {
  data.frame(gate = name, passed = passed, detail = detail, blocks = blocks,
             stringsAsFactors = FALSE)
}

#' Are the sex-chromosome genes even in the probe set?
#'
#' This has to be answered before the sex-gene gate means anything. A Flex
#' panel that does not target Xist or Ddx3y cannot detect them, which is a
#' documented limitation of the assay; the same absence when they ARE targeted
#' means the contrast is not measuring sex. The two look identical in a results
#' table and have opposite implications.
check_sex_genes_in_panel <- function(obj, cfg, assay = "RNA") {
  genes <- cfg$gates$sex_chromosome_genes
  present <- genes %in% rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  data.frame(gene = genes, in_panel = present, stringsAsFactors = FALSE)
}

#' Sex-chromosome genes must rank near the top of a male-vs-female contrast.
#'
#' Positive control. If Xist and the Y genes are not at the top, whatever the
#' contrast is measuring, it is not sex.
gate_sex_chromosome_genes <- function(ranked_genes, obj, cfg) {
  panel <- check_sex_genes_in_panel(obj, cfg)
  targeted <- panel$gene[panel$in_panel]

  if (!length(targeted))
    return(gate_result(
      "sex chromosome genes", NA,
      paste0("none of ", paste(panel$gene, collapse = ", "),
             " are in the probe set -- the panel cannot detect sex, which is a ",
             "limitation of the assay rather than a pipeline fault"),
      "any claim that the sex contrast is measuring sex"))

  top <- utils::head(ranked_genes, cfg$gates$sex_gene_expected_in_top)
  found <- intersect(targeted, top)
  ranks <- match(targeted, ranked_genes)

  gate_result(
    "sex chromosome genes", length(found) > 0,
    sprintf("%d of %d targeted sex genes in the top %d (ranks: %s)",
            length(found), length(targeted), cfg$gates$sex_gene_expected_in_top,
            paste(sprintf("%s=%s", targeted, ifelse(is.na(ranks), "absent", ranks)),
                  collapse = ", ")),
    "every cross-sex comparison")
}

#' Compared groups must be within a modest depth ratio.
#'
#' Evaluate this on the counts the analysis will actually use. Run against raw
#' counts it fails by construction -- the whole reason the matched assay exists
#' is that the raw libraries differ in depth -- so gating there would block the
#' pipeline on the problem the next step solves. `label` names which counts
#' were measured so the table says so.
gate_depth_ratio <- function(counts, groups, cfg, label = "compared groups") {
  ratio <- depth_ratio(counts, groups)
  gate_result(
    "depth ratio", ratio <= cfg$gates$max_depth_ratio,
    sprintf("%.2fx between %s (limit %.2fx); medians: %s",
            ratio, label, cfg$gates$max_depth_ratio,
            paste(sprintf("%s=%.0f", names(group_depth(counts, groups)),
                          group_depth(counts, groups)), collapse = ", ")),
    "any differential test on these counts")
}

#' No stratum below the minimum cell count may produce output.
#'
#' This gate's scope is the LISTED STRATA, not the step: a rare stage that no
#' run could ever populate is a reason to drop that stage, not to abandon the
#' other five. Callers that analyse strata one at a time should use
#' stratum_is_usable() inside the loop and reserve this gate for the question
#' it can actually answer -- whether anything is left to analyse at all.
gate_stratum_sizes <- function(counts_table, cfg, require_all = TRUE) {
  smallest <- min(counts_table)
  min_cells <- cfg$gates$min_cells_per_stratum
  failing <- names(which(counts_table < min_cells))
  usable <- names(which(counts_table >= min_cells))

  passed <- if (require_all) smallest >= min_cells else length(usable) > 0
  gate_result(
    "minimum cells per stratum", passed,
    sprintf("%d of %d strata at or above %d cells (smallest has %d)%s",
            length(usable), length(counts_table), min_cells, smallest,
            if (length(failing))
              paste0("; excluded: ", paste(utils::head(failing, 8), collapse = ", "))
            else ""),
    if (require_all) "stage-stratified results for the listed strata"
    else "every stage-stratified result, since no stratum is large enough")
}

#' Isotype controls must come out non-significant.
#'
#' Negative control on the ADT layer. If the isotypes separate the groups, the
#' protein layer is picking up a technical difference and cannot be trusted as
#' the depth-robust modality.
gate_isotype_controls <- function(isotype_results, cfg,
                                  effect_col = "effect_size", threshold = 0.25) {
  if (is.null(isotype_results) || !nrow(isotype_results))
    return(gate_result("ADT isotype controls", NA, "no isotype results supplied", ""))

  worst <- max(abs(isotype_results[[effect_col]]), na.rm = TRUE)
  gate_result(
    "ADT isotype controls", worst <= threshold,
    sprintf("largest isotype effect %.3f across %d controls (limit %.2f)",
            worst, nrow(isotype_results), threshold),
    "ADT-based stage assignment and any ADT contrast")
}

#' Run the gates, report them as a table, and stop if any fatal gate failed.
#'
#' NA means the gate could not be evaluated, which is reported but not fatal --
#' an unevaluated gate is not a passed gate, and the table says so.
run_gates <- function(gates, cfg, context = "") {
  results <- do.call(rbind, gates)
  results$status <- ifelse(is.na(results$passed), "NOT EVALUATED",
                    ifelse(results$passed, "pass", "FAIL"))

  log_step("control gates", if (nzchar(context)) paste0(" (", context, ")") else "", ":")
  for (i in seq_len(nrow(results)))
    log_step(sprintf("  [%-13s] %-26s %s", results$status[i], results$gate[i],
                     results$detail[i]))

  failed <- results[which(results$passed == FALSE), ]
  if (nrow(failed)) {
    message <- paste0(
      nrow(failed), " control gate(s) failed:\n",
      paste0("  - ", failed$gate, ": ", failed$detail,
             "\n    blocks: ", failed$blocks, collapse = "\n"),
      "\n\nThese are blocking because each marks a way this data has already been ",
      "observed to mislead.\nSet gates.fatal: false in config/config.yml to ",
      "continue anyway -- results produced that way are not reportable.")
    if (isTRUE(cfg$gates$fatal)) stop(message) else warning(message)
  }
  results
}
