#' Get time-related variable patterns
#'
#' @title Get time variable patterns
#'
#' @description Returns a regex pattern string containing common time-related variable names
#' used for detecting time variables in sample metadata.
#'
#' @return A character string containing pipe-separated time-related variable patterns
#'
#' @examples
#' get_timevars()
#' # Returns "second|minute|hour|day|week|month|year|time|age"
#'
#' @export
get_timevars <- function() {
  return(
    "second|minute|hour|day|week|month|year|time|age"
  )
}
