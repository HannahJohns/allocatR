#' Randomisation Schemes
#' @aliases RandomisationScheme
#' @description Functions that can create Randomisation Schemes to be used in
#' \link[=allocate]{\code{allocate()}} or \link[=trialSim]{\code{trialSim()}}
#' functions
#'
#' The following functions produce an object of type \code{RandomisationScheme}
#' @details
#' \itemize{
#' \item \link[=csmsb]{csmsb()}
#' }
#' @exportS3Method
print.RandomisationScheme <- function(x){

  print(x)

  return(invisible(x))
}
