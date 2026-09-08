##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Durable AI infographics
# =============================================================================
#
# Per-slot infographic images generated from an already-written AI report
# (pgx$ai[[slot]]$report), precomputed at compute time and stored so every
# viewer of a dataset sees the identical image instead of a per-session
# regeneration. Storage is a sibling of the report fields:
#
#   pgx$ai[[slot]]$infographic <- list(
#     status, error,           # "done" | "error"
#     bytes, content_type,     # raw image bytes -- the durable part
#     path,                    # per-session tempfile, session-local fallback
#     prompt, model,           # generation provenance
#     style, n_blocks, date, metadata
#   )
#
# v1 scope covers only the static report slots (combined, wgcna, wgcna_mox,
# mofa, de, pathways); dynamic drug-connectivity slots stay on-demand-only.

#' Read the pgx$ai slot as a plain list
#'
#' Shared by the report and infographic families (\code{ai_report_get()}/
#' \code{ai_report_slots()} on the omicsplayground side, \code{ai_infographic_get()}/
#' \code{ai_infographic_slots()} here) for a null-safe read of \code{pgx$ai}.
#'
#' @param pgx PGX object as a list.
#' @return \code{pgx$ai} as a list, or \code{NULL} when absent/malformed.
#' @export
ai_report_ai_slot <- function(pgx) {
  if (is.null(pgx) || !is.list(pgx)) return(NULL)
  ai <- pgx$ai
  if (is.null(ai) || !is.list(ai)) return(NULL)
  ai
}

#' Read a stored AI infographic for one report slot
#'
#' Returns an infographic entry only when it contains stored bytes, an existing
#' file path, or a stored error status.
#'
#' @param pgx PGX object as a list.
#' @param slot AI report slot name.
#' @return Infographic entry list, or \code{NULL} when unavailable.
#' @export
ai_infographic_get <- function(pgx, slot) {
  if (is.null(slot) || length(slot) != 1L) return(NULL)
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) NA_character_)
  if (is.na(slot) || !nzchar(slot)) return(NULL)

  ai <- ai_report_ai_slot(pgx)
  if (is.null(ai) || !is.list(ai[[slot]])) return(NULL)

  img <- ai[[slot]]$infographic
  if (!is.list(img)) return(NULL)
  has_bytes <- is.raw(img$bytes) && length(img$bytes) > 0L
  has_path <- is.character(img$path) && length(img$path) > 0L &&
    nzchar(img$path[[1L]]) && file.exists(img$path[[1L]])
  if (!has_bytes && !has_path && !identical(img$status, "error")) return(NULL)
  img
}

#' List AI report slots with stored infographic state
#'
#' @param pgx PGX object as a list.
#' @return Character vector of AI report slot names with infographic entries.
#' @export
ai_infographic_slots <- function(pgx) {
  ai <- ai_report_ai_slot(pgx)
  if (is.null(ai)) return(character(0))

  slots <- setdiff(names(ai), "meta")
  slots[vapply(slots, function(slot) {
    !is.null(ai_infographic_get(pgx, slot))
  }, logical(1))]
}

#' Convert provider errors into user-facing infographic messages
#'
#' @param error Error condition or character error message.
#' @return Sanitized message suitable for storage and UI display.
#' @export
ai_infographic_friendly_error <- function(error) {
  msg <- if (inherits(error, "condition")) {
    conditionMessage(error)
  } else if (is.null(error)) {
    ""
  } else {
    paste(as.character(error), collapse = " ")
  }
  msg <- trimws(msg)

  if (!nzchar(msg)) {
    return("Infographic generation failed. Please try again.")
  }
  if (grepl("503|Service Unavailable|high demand|overload|saturat",
      msg, ignore.case = TRUE)) {
    return(paste0(
      "Image generation server is temporarily overloaded. ",
      "Please try again in a few minutes."
    ))
  }
  if (grepl("429|Too Many Requests|rate.?limit", msg, ignore.case = TRUE)) {
    return("Image generation rate limit reached. Please wait a moment and try again.")
  }
  if (grepl("timeout|timed.?out", msg, ignore.case = TRUE)) {
    return(paste0(
      "Image generation timed out. The service may be busy; ",
      "please try again."
    ))
  }
  if (grepl("No image data|Empty response|no response", msg, ignore.case = TRUE)) {
    return(paste0(
      "Image generation server seems saturated and returned no image. ",
      "Please try again."
    ))
  }
  if (grepl("_API_KEY|API.?key|401|Unauthorized", msg, ignore.case = TRUE)) {
    return(paste0(
      "Image generation is not configured correctly. ",
      "Please contact your administrator."
    ))
  }

  "Image generation failed. Please try again."
}

#' Extract persisted image payload fields from an omicsai image result
#'
#' @param result \code{omicsai_image_result} list, or \code{NULL} for errors.
#' @return List containing bytes, content type, path, prompt, model, and metadata.
#' @export
ai_infographic_payload <- function(result) {
  path <- NULL
  prompt <- NULL
  metadata <- NULL
  if (is.list(result)) {
    path <- if (is.character(result$path) && length(result$path) > 0L) {
      result$path[[1L]]
    } else {
      NULL
    }
    prompt <- result$prompt
    metadata <- result$metadata
  }

  bytes <- raw(0)
  if (!is.null(path) && file.exists(path)) {
    bytes <- readBin(path, what = raw(), n = file.info(path)$size)
  }

  ext <- tolower(tools::file_ext(if (is.null(path)) "image.png" else path))
  list(
    bytes = bytes,
    content_type = switch(ext,
      jpg = "image/jpeg",
      jpeg = "image/jpeg",
      "image/png"),
    path = path,
    prompt = prompt,
    model = if (is.null(metadata$model)) NULL else metadata$model,
    metadata = metadata
  )
}

#' Store one AI infographic under pgx$ai
#'
#' Writes generated image bytes and metadata, or a sanitized error entry, into
#' \code{pgx$ai[[slot]]$infographic}.
#'
#' @param pgx PGX object as a list.
#' @param slot AI report slot name.
#' @param result \code{omicsai_image_result} list, or \code{NULL} for errors.
#' @param status Infographic status, usually \code{"done"} or \code{"error"}.
#' @param error Optional provider error for failed generation.
#' @param style omicsai image style ID.
#' @param n_blocks Deprecated layout metadata retained for compatibility.
#' @return Updated PGX object.
#' @export
ai_infographic_set <- function(pgx, slot, result,
                               status = "done", error = NULL,
                               style = NULL, n_blocks = NULL) {
  if (is.null(pgx) || !is.list(pgx)) return(pgx)
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) NA_character_)
  if (is.na(slot) || !nzchar(slot)) return(pgx)

  if (is.null(pgx$ai) || !is.list(pgx$ai)) pgx$ai <- list()
  if (is.null(pgx$ai[[slot]]) || !is.list(pgx$ai[[slot]])) {
    pgx$ai[[slot]] <- list()
  }
  # Normalize generated image data before writing the PGX AI schema entry.
  payload <- ai_infographic_payload(result)

  pgx$ai[[slot]]$infographic <- list(
    status = status,
    error = if (identical(status, "error")) {
      ai_infographic_friendly_error(error)
    } else {
      error
    },
    bytes = payload$bytes,
    content_type = payload$content_type,
    path = payload$path,
    prompt = payload$prompt,
    model = payload$model,
    style = style,
    n_blocks = n_blocks,
    date = as.character(Sys.time()),
    metadata = payload$metadata
  )
  pgx
}

# -----------------------------------------------------------------------------
# Prompt assembly
# -----------------------------------------------------------------------------

# Static v1 label map for the report slots eligible for precomputed
# infographics. Drug-connectivity slots (drugs_*) are app-layer and excluded.
.AI_INFOGRAPHIC_LABELS <- c(
  combined = "Summary",
  wgcna = "WGCNA",
  wgcna_mox = "moxWGCNA",
  mofa = "MOFA",
  de = "Differential Expression",
  pathways = "Enrichment"
)

# Build the system + board omicsai prompt for one slot's infographic. Mirrors
# infographic_build_image_prompt() (infographic_server.R), with label
# resolution simplified to the static v1 slot map (no drug-slot branch).
.ai_infographic_build_prompt <- function(report, slot, pgx, style) {
  label <- .AI_INFOGRAPHIC_LABELS[[slot]] %||% slot
  organism <- pgx$organism
  if (is.null(organism) || !nzchar(organism)) organism <- "human"
  clean <- omicsai::omicsai_strip_report_noise(report)
  species <- omicsai::omicsai_image_species_visual(organism)
  prompt <- omicsai::image_prompt(
    role = omicsai::frag("system_base"),
    task = omicsai::frag("image/infographic",
      params = list(board_name = label)),
    species = species,
    style = omicsai::frag(paste0("image/styles/", style)),
    report = clean
  )
  omicsai::build_prompt(prompt)
}

# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

# v1 static slot list eligible for precomputed infographics.
.AI_INFOGRAPHIC_SLOTS <- c("combined", "wgcna", "wgcna_mox", "mofa", "de", "pathways")

#' Precompute durable AI infographics for a PGX object.
#'
#' Iterates the v1 static report slots and, for each slot with an existing
#' \code{pgx$ai[[slot]]$report}, generates an infographic image via
#' \code{omicsai::omicsai_gen_image()} and stores it under
#' \code{pgx$ai[[slot]]$infographic}. Slots without a report, or slots that
#' already carry an infographic (unless \code{ai$force}), are skipped.
#'
#' @param pgx PGX object.
#' @param ai Infographic-generation options; \code{NULL} returns \code{pgx}
#'   unchanged. Recognised fields: \code{style}, \code{img_model},
#'   \code{select} (slot names, defaults to the v1 static list),
#'   \code{force}, \code{on_error}.
#' @return PGX object with \code{pgx$ai[[slot]]$infographic} populated.
#' @export
pgx.update_infographics <- function(pgx, ai = NULL) {
  if (is.null(ai)) return(pgx)
  ai <- .ai_resolve_defaults(ai)
  if (is.null(pgx$ai) || !is.list(pgx$ai)) pgx$ai <- list()

  select <- ai$select %||% .AI_INFOGRAPHIC_SLOTS
  select <- intersect(select, .AI_INFOGRAPHIC_SLOTS)
  style <- ai$style %||% "bigomics"

  for (slot in select) {
    if (is.null(pgx$ai[[slot]]$report)) {
      info("[pgx.update_infographics] '", slot, "' has no report; skipping")
      next
    }
    if (!is.null(pgx$ai[[slot]]$infographic) && !isTRUE(ai$force)) {
      info("[pgx.update_infographics] '", slot, "' already present; skipping")
      next
    }

    info("[pgx.update_infographics] generating '", slot, "' infographic...")
    pgx <- tryCatch(
      {
        built <- .ai_infographic_build_prompt(pgx$ai[[slot]]$report, slot, pgx, style)
        config <- omicsai::omicsai_image_config(
          model = ai$img_model,
          system_prompt = built$system,
          style = style,
          image_size = "1K"
        )
        result <- omicsai::omicsai_gen_image(
          template = built$board,
          params = NULL,
          config = config
        )
        ai_infographic_set(pgx, slot, result, status = "done", style = style)
      },
      error = function(e) {
        msg <- paste0("[pgx.update_infographics] '", slot, "' failed: ",
                      conditionMessage(e))
        if (identical(ai$on_error, "abort")) stop(msg, call. = FALSE)
        if (identical(ai$on_error, "warn"))  warning(msg, call. = FALSE)
        else                                  info(msg)
        ai_infographic_set(pgx, slot, result = NULL, status = "error",
          error = e, style = style)
      }
    )
  }

  pgx
}
