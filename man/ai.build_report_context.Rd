\name{ai.build_report_context}
\alias{ai.build_report_context}
\title{Build compact text context for Copilot}
\usage{
ai.build_report_context(
  pgx,
  sections = NULL,
  collate = TRUE,
  ntop = 20,
  max_chars = Inf,
  ...
)
}
\arguments{
\item{pgx}{PGX object or list.}

\item{sections}{Optional character vector of context sections to include.}

\item{collate}{Logical; return one text string when \code{TRUE}, or a named list when \code{FALSE}.}

\item{ntop}{Number of rows/items to show for tabular or list sections.}

\item{max_chars}{Maximum number of characters to return. Defaults to \code{Inf}.}

\item{...}{Reserved for compatibility.}
}
\value{
A character string when \code{collate = TRUE}, otherwise a named list of character strings.
}
\description{
Builds Copilot context from dataset metadata and precomputed reports stored in
\code{pgx$ai}.
}
