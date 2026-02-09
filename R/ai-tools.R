
##======================================================================
##==================== TOOLS ===========================================
##======================================================================

#' Plot Volcano
#'
#' @export
ai.tool_plot_volcano <- ellmer::tool(
  function(contrast, psig = 0.05) {
    message("calling tool: plot_volcano(). contrast = ", contrast)
    cmd <- "playbase::pgx.Volcano(pgx, contrast=1, psig=0.05)"
    paste0("<code>",cmd,"</code>")
  },
  "Code using playbase library to create volcano plot for a given contrast",
  contrast = ellmer::type_string(
    "The comparison/contrast for the Volcano plot",
    required = TRUE
  ),
  psig = ellmer::type_number(
    description = "Significance level. Default psig=0.05",
    required = FALSE
  )
)

#' Get expression values for gene
#'
#' @export
ai.tool_get_expression <- ellmer::tool(
  function(gene) {
    message("calling tool: get_expression(). gene = ", gene)
    pgx <- playdata::GEIGER_PGX
    list(
      values = pgx$X[gene,],
      plot_command = "<code>base::barplot(pgx$X[gene,])</code>"
    )
  },
  "Get expression values for a gene",
  gene = ellmer::type_string(
    description = "The gene name to retrieve the expression values for",
    required = TRUE
  )
)

#' Gets the current time in the given time zone.
#'
#' @export
ai.tool_get_current_time <- ellmer::tool(
  function(tz = "UTC") {
    message("calling tool: get_current_time()")
    format(Sys.time(), tz = tz, usetz = TRUE)
  },
  "Gets the current time in the given time zone.",
  tz = ellmer::type_string(
    description = "The time zone to get the current time in. Defaults to `\"UTC\"`.",
    required = FALSE
  )
)

