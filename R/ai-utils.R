##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Export markdown to PDF
#'
#' This keeps the historical Quarto-based report export flow used by the
#' webapp, while moving argument and output checks into small private helpers.
#' The intent is to preserve the familiar rendering logic and make failures
#' easier to diagnose when report export runs in production.
#'
#' @param text Character vector containing markdown content to render. Multiple
#'   values are collapsed with newlines before rendering.
#' @param file Single output path ending in `.pdf` or `.docx`.
#' @param tmpdir Optional temporary directory used for the intermediate Quarto
#'   project files. Defaults to `tempdir()`.
#' @param engine Quarto/LaTeX PDF engine name, such as `"pdflatex"` or
#'   `"lualatex"`.
#' @param font Optional font or LaTeX font family passed into Quarto frontmatter.
#'   When `NULL`, a default is chosen for supported PDF engines.
#' @param documentclass LaTeX document class to use for PDF output.
#' @param force.ascii Logical flag indicating whether markdown content should be
#'   converted to ASCII before rendering.
#' @param logo Optional path to a logo image added to the PDF header when the
#'   file exists.
#' @param quiet Logical flag passed to `quarto::quarto_render()` to suppress
#'   Quarto output.
#' @export
markdownToPDF <- function(text, file, tmpdir = NULL, engine = "pdflatex",
                          font = "mathpazo", documentclass = "report",
                          force.ascii = TRUE, logo = NULL, quiet = TRUE) {
  # Normalize inputs and collect render paths before building Quarto frontmatter.
  x <- .markdownToPDF_prepare(
    text, file, tmpdir, engine, documentclass, force.ascii, logo
  )

  text <- x$text
  file <- x$file
  tmpdir <- x$tmpdir
  logo <- x$logo
  ext <- x$ext

  if (is.null(font) && engine == "pdflatex") font <- "mathpazo"
  if (is.null(font) && engine == "lualatex") font <- "Lato"

  if (engine == "lualatex") {
    text <- gsub("[.]underline", "", text)
  }

  ## header/frontmatter
  hdr <- c(
    "---",
    "format:",
    "  pdf:",
    paste0("    pdf-engine: ", engine),
    paste0("    documentclass: ", documentclass),
    "    papersize: a4",
    "    fontsize: 10pt",
    "    geometry:",
    "      - left=24mm",
    "      - right=20mm",
    "      - top=25mm",
    "      - bottom=18mm"
  )
  if (!is.null(font) && engine == "lualatex") {
    hdr <- c(hdr, paste0("    mainfont: ", font))
  }
  if (!is.null(font) && engine == "pdflatex") {
    hdr <- c(hdr, paste0("    fontfamily: ", font))
  }
  hdr <- c(hdr, "    fig-pos: 'h!'")
  if (!is.null(logo)) {
    tex_header <- tempfile(fileext = ".tex")
    writeLines(c(
      "\\usepackage{graphicx}",
      "\\usepackage{eso-pic}",
      paste0(
        "\\AddToShipoutPictureBG*{\\put(\\LenToUnit{\\paperwidth-20mm},",
        "\\LenToUnit{\\paperheight-24mm}){\\makebox[0pt][r]{",
        "\\includegraphics[height=9mm]{",
        normalizePath(logo, winslash = "/", mustWork = TRUE),
        "}}}}"
      )
    ), tex_header)
    hdr <- c(hdr, paste0(
      "    include-in-header: ",
      shQuote(normalizePath(tex_header, winslash = "/", mustWork = TRUE))
    ))
  }
  hdr <- paste0(paste(c(hdr, "---"), collapse = "\n"), "\n\n")
  text <- paste0(hdr, text)

  workdir <- tempfile("markdownToPDF-", tmpdir = tmpdir)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(workdir)) {
    stop("Unable to create temporary render directory: ", workdir, call. = FALSE)
  }
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  qmd.file <- file.path(workdir, "report.qmd")
  writeLines(text, con = qmd.file, useBytes = TRUE)

  out.file <- file.path(workdir, paste0("report.", ext))
  quarto::quarto_render(
    qmd.file,
    output_format = ext,
    output_file = basename(out.file),
    quarto_args = c("--output-dir", workdir),
    quiet = quiet
  )

  # Validate the rendered artifact and move it to the requested output path.
  .markdownToPDF_copy(out.file, file)
}

# Keep input validation and path setup outside the render body so failures are
# reported before Quarto starts and the main function stays focused on rendering.
.markdownToPDF_prepare <- function(text, file, tmpdir, engine, documentclass,
                                   force.ascii, logo) {
  if (missing(text) || is.null(text)) {
    stop("'text' must contain markdown content.", call. = FALSE)
  }
  if (missing(file) || length(file) != 1 || is.na(file) || file == "") {
    stop("'file' must be a single output path.", call. = FALSE)
  }
  if (length(engine) != 1 || is.na(engine) || engine == "") {
    stop("'engine' must be a single Quarto PDF engine name.", call. = FALSE)
  }
  if (length(documentclass) != 1 || is.na(documentclass) || documentclass == "") {
    stop("'documentclass' must be a single LaTeX document class.", call. = FALSE)
  }
  text <- paste(as.character(text), collapse = "\n")
  text <- gsub(intToUtf8(8209), "-", text, fixed = TRUE)
  if (isTRUE(force.ascii)) text <- iconv2ascii(text)
  text <- gsub("```markdown|```", "", text)

  file <- path.expand(file)
  ext <- tolower(tools::file_ext(file))
  if (!ext %in% c("pdf", "docx")) {
    stop("'file' must end in '.pdf' or '.docx'.", call. = FALSE)
  }
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("[markdownToPDF] ERROR: package 'quarto' is required to export markdown reports",
         call. = FALSE)
  }

  outdir <- dirname(file)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(outdir)) {
    stop("Unable to create output directory: ", outdir, call. = FALSE)
  }

  if (is.null(tmpdir)) tmpdir <- tempdir()
  tmpdir <- path.expand(tmpdir)
  if (!dir.exists(tmpdir)) {
    dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(tmpdir)) {
    stop("Unable to create temporary directory: ", tmpdir, call. = FALSE)
  }

  if (!is.null(logo)) {
    logo <- path.expand(logo)
    if (!file.exists(logo)) logo <- NULL
  }

  list(text = text, file = file, tmpdir = tmpdir, logo = logo, ext = ext)
}

# Keep final artifact checks in one helper so every render verifies that Quarto
# produced a non-empty file before replacing the caller's requested destination.
.markdownToPDF_copy <- function(out.file, file) {
  out.info <- file.info(out.file)
  if (!file.exists(out.file) || is.na(out.info$size) || out.info$size == 0) {
    stop("Quarto did not create a non-empty report file: ", out.file,
         call. = FALSE)
  }
  if (!isTRUE(file.copy(out.file, file, overwrite = TRUE))) {
    stop("Unable to copy rendered report to: ", file, call. = FALSE)
  }
  file
}
