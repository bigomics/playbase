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
  tier         = NULL,
  credentials  = NULL,
  ## A report module routinely decodes for 120-170s, so the per-request
  ## deadline has to clear that comfortably; two attempts caps the worst case
  ## at 480s rather than ellmer's default 900s.
  timeout_seconds = 240L,
  retries         = 2L,
  ## Reports are a summarisation job over data that has already been computed,
  ## so extended reasoning buys little and costs a lot of wall time: on a
  ## reasoning model at "medium" the reasoning tokens were ~60% of the latency
  ## for a report no longer than the one produced without them. Models that do
  ## not accept the knob ignore it.
  reasoning_effort = "low"
)

.ai_resolve_defaults <- function(ai) {
  if (is.null(ai)) ai <- list()
  if (!is.list(ai)) stop("[.ai_resolve_defaults] `ai` must be a list")
  out <- .AI_DEFAULTS
  for (nm in names(ai)) out[[nm]] <- ai[[nm]]
  if (!out$report_type %in% c("normal", "deep")) {
    stop("[.ai_resolve_defaults] ai$report_type must be 'normal' or 'deep'; got '",
         out$report_type, "'")
  }
  if (!out$on_error %in% c("skip", "abort", "warn")) {
    stop("[.ai_resolve_defaults] ai$on_error must be 'skip', 'abort', or 'warn'")
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

.ai_report_normheadings <- function(report) {
  omicsai::omicsai_normheadings(report)
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

.ai_report_module_builder <- function(module) {
  fn <- .ai_report_module_function(module)
  if (is.null(fn)) NULL else sub("\\.create_report$", ".build_jobs", fn)
}

#' Describe one LLM call for a report slot
#'
#' A job is the unit the orchestrator hands to a runner. It deliberately holds
#' nothing but plain strings plus a `finalize` closure, so the two halves can be
#' separated: `system`/`board` are all a worker process needs (a few hundred KB,
#' versus the hundreds of MB a pgx costs to ship), while `finalize` stays with
#' the caller and never crosses a process boundary.
#'
#' @param module Module name the job belongs to.
#' @param slot Target slot under `pgx$ai`.
#' @param bp Built prompt, `list(system=, board=)`.
#' @param finalize Function applied to the generated report text, or NULL.
#' @return An `ai_report_job` list.
#' @keywords internal
.ai_report_job <- function(module, slot, bp, finalize = NULL) {
  structure(
    list(
      module   = module,
      slot     = slot,
      system   = bp$system,
      board    = bp$board,
      finalize = finalize
    ),
    class = c("ai_report_job", "list")
  )
}

#' Build the LLM jobs for one module
#'
#' @return A list of `ai_report_job`s, or an `ai_report_skip`.
#' @keywords internal
.ai_report_build_module <- function(module, pgx, ai) {
  fn_name <- .ai_report_module_builder(module)
  if (is.null(fn_name)) {
    return(structure(list(reason = "module not found"),
                     class = "ai_report_skip"))
  }
  if (!exists(fn_name, mode = "function")) {
    return(structure(list(reason = paste0("no entry point '", fn_name, "()'")),
                     class = "ai_report_skip"))
  }
  slice <- .ai_report_module_slice(module, pgx)
  if (is.null(slice)) {
    return(structure(list(reason = "slot is empty"), class = "ai_report_skip"))
  }
  get(fn_name, mode = "function")(pgx, slice, ai)
}

#' Run one report job against the provider
#'
#' @param job An `ai_report_job`.
#' @param ai Resolved AI options.
#' @return `list(report=, prompt=, usage=)`.
#' @keywords internal
.ai_report_run_job <- function(job, ai) {
  out <- .ai_report_run_prompt(list(system = job$system, board = job$board), ai)
  if (is.function(job$finalize)) out$report <- job$finalize(out$report)
  out
}

#' Legacy single-module entry point, expressed as build + run
#'
#' Kept so `ai.<module>.create_report()` callers outside the orchestrator keep
#' working unchanged.
#' @keywords internal
.ai_report_create_report_compat <- function(module, pgx, slice, ai) {
  ## Re-dispatch through the builder rather than the slice we were handed: the
  ## builders read the slice back off pgx themselves, and wgcna_mox may
  ## redirect to the single-omics builder.
  jobs <- .ai_report_build_module(module, pgx, ai)
  if (inherits(jobs, "ai_report_skip") || !length(jobs)) return(NULL)
  .ai_report_run_job(jobs[[1L]], ai)
}

#' Build every LLM job needed for the selected report modules
#'
#' Assembles prompts without contacting any provider, so a caller can run them
#' how it likes - serially, or fanned out across worker processes. Prompt
#' assembly is cheap (single-digit seconds for a full pgx); essentially all the
#' wall time of a report run is the provider calls these jobs stand for.
#'
#' `combined` is deliberately *not* returned here even when selected: its
#' prompt is assembled from `pgx$ai`, so it can only be built once the other
#' modules' reports have been folded back in. Call
#' [pgx.build_combined_report_job()] for that second phase.
#'
#' @param pgx PGX object.
#' @param ai Report-generation options, as for [pgx.update_reports()].
#' @return A list of `ai_report_job`s (possibly empty).
#' @export
pgx.build_report_jobs <- function(pgx, ai = NULL) {
  ai <- .ai_resolve_defaults(ai)
  modules <- setdiff(ai$select, "combined")
  jobs <- list()
  for (module in modules) {
    if (!.ai_report_wanted(pgx, module, ai)) next
    built <- tryCatch(
      .ai_report_build_module(module, pgx, ai),
      error = function(e) {
        .ai_report_report_error(module, e, ai)
        NULL
      }
    )
    if (inherits(built, "ai_report_skip")) {
      message("[pgx.build_report_jobs] skipping '", module, "', ", built$reason)
      next
    }
    if (!length(built)) next
    jobs <- c(jobs, built)
  }
  jobs
}

#' Build the combined-summary job for the current state of `pgx$ai`
#'
#' @inheritParams pgx.build_report_jobs
#' @return A list holding zero or one `ai_report_job`.
#' @export
pgx.build_combined_report_job <- function(pgx, ai = NULL) {
  ai <- .ai_resolve_defaults(ai)
  if (!"combined" %in% ai$select) return(list())
  if (!.ai_report_wanted(pgx, "combined", ai)) return(list())
  built <- tryCatch(
    .ai_report_build_module("combined", pgx, ai),
    error = function(e) {
      .ai_report_report_error("combined", e, ai)
      NULL
    }
  )
  if (inherits(built, "ai_report_skip")) {
    message("[pgx.build_combined_report_job] skipping combined, ", built$reason)
    return(list())
  }
  if (!length(built)) list() else built
}

#' Fold one finished report into `pgx$ai`
#'
#' @param pgx PGX object.
#' @param job The `ai_report_job` the result came from.
#' @param result `list(report=, prompt=, usage=)` as produced by a runner. The
#'   `report` is passed through the job's `finalize` unless `finalized = TRUE`.
#' @param finalized Set TRUE when the runner already applied `job$finalize`.
#' @return The PGX object with the slot populated.
#' @export
pgx.apply_report_result <- function(pgx, job, result, finalized = FALSE) {
  if (is.null(result) || is.null(result$report)) return(pgx)
  report <- result$report
  if (!isTRUE(finalized) && is.function(job$finalize)) {
    report <- job$finalize(report)
  }
  if (is.null(pgx$ai)) pgx$ai <- list()
  pgx$ai[[job$slot]] <- list(
    report     = .ai_report_normheadings(report),
    prompt     = result$prompt,
    usage      = result$usage,
    created_at = as.numeric(Sys.time()),
    edited     = FALSE,
    edited_at  = ""
  )
  pgx
}

#' Should this module be (re)generated?
#' @keywords internal
.ai_report_wanted <- function(pgx, module, ai) {
  if (isTRUE(ai$force)) return(TRUE)
  if (is.null(pgx$ai[[module]])) return(TRUE)
  message("[pgx.update_reports] '", module, "' already present; skipping")
  FALSE
}

#' Apply the configured `on_error` policy to a module failure.
#' @keywords internal
.ai_report_report_error <- function(module, e, ai) {
  msg <- paste0("[pgx.update_reports] '", module, "' failed: ",
                conditionMessage(e))
  if (identical(ai$on_error, "abort")) stop(msg, call. = FALSE)
  if (identical(ai$on_error, "warn"))  warning(msg, call. = FALSE)
  else                                  message(msg)
  NULL
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

## Superseded by .ai_report_build_module() + .ai_report_run_job(), which split
## prompt assembly from the provider call so the two can run in different
## processes. Retained because it is the documented extension point for a
## module entry point with the (pgx, slice, ai) signature.
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

  ## Two phases, because `combined` reads pgx$ai: everything else is built and
  ## run first, then combined is built against the updated pgx. Within a phase
  ## the jobs are run serially here; callers that want them concurrent build
  ## the same jobs with pgx.build_report_jobs() and run them themselves.
  run_phase <- function(pgx, jobs) {
    for (job in jobs) {
      message("[pgx.update_reports] generating '", job$slot, "' report...")
      res <- tryCatch(
        .ai_report_run_job(job, ai),
        error = function(e) .ai_report_report_error(job$module, e, ai)
      )
      if (is.null(res)) {
        message("[pgx.update_reports] skipping '", job$slot,
                "', no report returned")
        next
      }
      pgx <- pgx.apply_report_result(pgx, job, res, finalized = TRUE)
    }
    pgx
  }

  pgx <- run_phase(pgx, pgx.build_report_jobs(pgx, ai))
  pgx <- run_phase(pgx, pgx.build_combined_report_job(pgx, ai))

  pgx$ai$meta <- .ai_build_meta(ai)
  pgx
}
