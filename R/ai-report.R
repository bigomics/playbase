##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# AI report orchestrator
# =============================================================================

.AI_DEFAULTS <- list(
  llm_model    = "openai:gpt-5.4-mini",
  img_model    = NULL,
  report_type  = "normal",
  select       = c("wgcna", "wgcna_mox", "mofa", "drugs", "de", "pathways",
                   "combined"),
  drug_dbs     = NULL,
  ntop         = 100,
  psig         = 0.05,
  userprompt   = NULL,
  force        = FALSE,
  on_error     = "skip",
  max_turns    = 50L,
  tier         = NULL
)

.ai_resolve_defaults <- function(ai) {
  if (is.null(ai)) ai <- list()
  if (!is.list(ai)) stop("[pgx.update_reports] `ai` must be a list")
  out <- .AI_DEFAULTS
  for (nm in names(ai)) out[[nm]] <- ai[[nm]]
  if (!out$report_type %in% c("normal", "deep")) {
    stop("[pgx.update_reports] ai$report_type must be 'normal' or 'deep'; got '",
         out$report_type, "'")
  }
  if (!out$on_error %in% c("skip", "abort", "warn")) {
    stop("[pgx.update_reports] ai$on_error must be 'skip', 'abort', or 'warn'")
  }
  out
}

.ai_build_meta <- function(ai) {
  omicsai_version <- tryCatch(
    as.character(utils::packageVersion("omicsai")),
    error = function(e) NA_character_
  )
  list(
    date            = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    llm_model       = ai$llm_model,
    img_model       = ai$img_model,
    report_type     = ai$report_type,
    select          = ai$select,
    drug_dbs        = ai$drug_dbs,
    ntop            = ai$ntop,
    psig            = ai$psig,
    omicsai_version = omicsai_version,
    options         = list()
  )
}

.ai_report_module_function <- function(module) {
  switch(module,
    wgcna     = "ai.wgcna.create_report",
    wgcna_mox = "ai.wgcna_mox.create_report",
    mofa      = "ai.mofa.create_report",
    drugs     = "ai.drugs.create_report",
    de        = "ai.de.create_report",
    pathways  = "ai.pathways.create_report",
    combined  = "ai.combined.create_report",
    NULL
  )
}

.ai_report_module_slice <- function(module, pgx) {
  switch(module,
    wgcna     = pgx$wgcna,
    wgcna_mox = pgx$wgcna_mox,
    mofa      = pgx$mofa,
    drugs     = pgx$drugs,
    de        = pgx$gx.meta,
    pathways  = pgx$gset.meta,
    combined  = pgx$ai,
    NULL
  )
}

.ai_dispatch_module <- function(module, pgx, ai) {
  fn_name <- .ai_report_module_function(module)
  if (is.null(fn_name)) {
    return(structure(list(reason = "module not found"),
                     class = "ai_report_skip"))
  }
  if (!exists(fn_name, mode = "function")) {
    return(structure(list(reason = paste0("no entry point '", fn_name, "()'")),
                     class = "ai_report_skip"))
  }
  fn <- get(fn_name, mode = "function")
  fn_args <- names(formals(fn))
  if (length(fn_args) < 3 || !identical(fn_args[1:3], c("pgx", "slice", "ai"))) {
    return(structure(list(
      reason = paste0("'", fn_name, "' does not use signature (pgx, slice, ai)")
    ), class = "ai_report_skip"))
  }
  slice <- .ai_report_module_slice(module, pgx)
  if (is.null(slice)) {
    return(structure(list(reason = "slot is empty"), class = "ai_report_skip"))
  }
  fn(pgx, slice, ai)
}

#' Generate AI reports for a PGX object.
#'
#' Dispatches each selected module to `ai.<module>.create_report(pgx, slice, ai)`
#' and stores reports under `pgx$ai$<module>$report` and
#' `pgx$ai$<module>$prompt`. Run-level metadata is stored in `pgx$ai$meta`.
#'
#' @param pgx PGX object.
#' @param ai Report-generation options. `NULL` returns `pgx` unchanged.
#' @return PGX object with `pgx$ai` populated.
#' @export
pgx.update_reports <- function(pgx, ai = NULL) {
  if (is.null(ai)) return(pgx)
  ai <- .ai_resolve_defaults(ai)

  if (identical(ai$report_type, "deep")) {
    stop("[pgx.update_reports] ai$report_type = 'deep' is not implemented. ",
         "Use ai$report_type = 'normal' for now.")
  }

  force_combined_only <- isTRUE(ai$force) &&
    length(ai$select) == 1L &&
    identical(ai$select[[1]], "combined")
  if (force_combined_only && !is.null(pgx$ai)) {
    pgx$ai$combined <- NULL
    ai$force <- FALSE
  } else if (isTRUE(ai$force)) {
    pgx$ai <- NULL
  }
  if (is.null(pgx$ai)) pgx$ai <- list()

  for (module in ai$select) {
    if (!is.null(pgx$ai[[module]]) && !isTRUE(ai$force)) {
      message("[pgx.update_reports] '", module, "' already present; skipping")
      next
    }
    message("[pgx.update_reports] generating '", module, "' report...")
    out <- tryCatch(
      .ai_dispatch_module(module, pgx, ai),
      error = function(e) {
        msg <- paste0("[pgx.update_reports] '", module, "' failed: ",
                      conditionMessage(e))
        if (identical(ai$on_error, "abort")) stop(msg, call. = FALSE)
        if (identical(ai$on_error, "warn"))  warning(msg, call. = FALSE)
        else                                  message(msg)
        NULL
      }
    )
    if (inherits(out, "ai_report_skip")) {
      message("[pgx.update_reports] skipping '", module, "', ", out$reason)
      next
    }
    if (is.null(out)) {
      message("[pgx.update_reports] skipping '", module,
              "', no report returned")
      next
    }
    if (inherits(out, "ai_report_multi")) {
      for (key in names(out)) {
        slot_name <- paste0(module, "_", key)
        pgx$ai[[slot_name]] <- list(report = out[[key]]$report,
                                    prompt = out[[key]]$prompt)
      }
    } else {
      pgx$ai[[module]] <- list(report = out$report, prompt = out$prompt)
    }
  }

  pgx$ai$meta <- .ai_build_meta(ai)
  pgx
}
