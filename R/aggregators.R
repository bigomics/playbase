#' Plot data aggregators
#'
#' The data aggregators are usefult to downsample data before plotting. This becomes
#' very helpful when plotting large datasets. The aggregators are used by the 
#' shinyHugePlot::downsample objects. Currently, the following aggregators are
#' random aggregator and gaussian aggregator. The random aggregator select the 
#' data following a uniform probability distribution. The gaussian aggregator
#' assumes the data is a multivariate gaussian distribution and select the data
#' from the less dense areas in order to preserve the information.
#' 
#' @param x Numeric vector of x data
#' @param y Numeric vector of y data  
#' @param n_out Integer specifying number of data points to sample
#' @docType class
#' @format An \code{R6::R6Class} object
#' @export
random_aggregator <- R6::R6Class(
  "random_aggregator",
  inherit = shinyHugePlot::aggregator,
  public = list(
    #' @description
    #' Constructor of the Aggregator.
    #' @param interleave_gaps,coef_gap,NA_position,accepted_datatype,...
    #' Arguments pass to the constructor of \code{aggregator} object.
    initialize = function(
      ...,
      interleave_gaps, coef_gap, NA_position, accepted_datatype
    ) {
      args <- c(as.list(environment()), list(...))
      do.call(super$initialize, args)
    }
  ),
  private = list(
    aggregate_exec = function(x, y, n_out) {

      idx <- sample(seq_along(x), size = n_out, replace = FALSE)

      return(list(x = x[idx], y = y[idx]))

    }
  )
)

#' @docType class
#' @format An \code{R6::R6Class} object
#' @describeIn random_aggregator gaussian aggregator constructor for the downsampler object
#' @export
gaussian_aggregator <- R6::R6Class(
  "gaussian_aggregator",
  inherit = shinyHugePlot::aggregator,
  public = list(
    #' @description
    #' Constructor of the Aggregator.
    #' @param interleave_gaps,coef_gap,NA_position,accepted_datatype,...
    #' Arguments pass to the constructor of \code{aggregator} object.
    initialize = function(
      ...,
      interleave_gaps, coef_gap, NA_position, accepted_datatype
    ) {
      args <- c(as.list(environment()), list(...))
      do.call(super$initialize, args)
    }
  ),
  private = list(
    aggregate_exec = function(x, y, n_out) {

      m <- matrix(c(x, y), ncol = 2)
      mus <- apply(m, 2, mean)
      sigmas <- cov(m)
      dens <- mvtnorm::dmvnorm(m, mean = mus, sigma = sigmas)
      sampe_p <- 1/dens
      idx <- sample(seq_along(x), size = n_out, replace = FALSE,  prob = sampe_p)

      return(list(x = x[idx], y = y[idx]))

    }
  )
)
