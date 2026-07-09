#' Randomise a participant
#'
#' @description Randomise a participant according to the MSB, CS-MSB or CSSize-MSB algorithms.
#'
#' @usage allocate( method = c("simple","csmsb","msb"),
#'                  arm,
#'                  data = NULL,
#'                  newData = NULL,
#'                  settings = NULL
#'                  allocationSequence = NULL,
#'                  implementation = c("cpp", "R")
#'                  )
#'
#' @param method A string indicating the randomisation method to be used
#' @param arm An integer vector containing the current treatment assignments
#' @param data A data frame with \code{length(arm)} rows containing relevant data (e.g. covariates) for currently assigned participants
#' @param newData A one-row data frame with the same columns as \code{data} containning relevant data (e.g. covariates) for the new participant to be assigned
#' @param settings A list containing settings to be used by a specific randomisation method
#' @param allocationSequence An object of type \code{Alloc_seq}
#' @param implementation A string indicating which implementation of a given method should be used
#'
#' @returns A list with class \code{Allocation} containing the following:
#' \describe{
#' \item assignment
#' \item details a list containing internal steps taken to reach this assignment result
#' }
#'
#' @details
#' Additional details...
#'
#'
allocate <- function(  method = c("simple","csmsb","msb"),
                       arm,
                       randSequence,
                       data = NULL,
                       newData = NULL,
                       settings = NULL,
                       implementation = NULL
){

  # Generic input sanitation
  data <- as.data.frame(data)
  newData <- as.data.frame(newData)

  if(method == "csmsb"){

    if(is.null(implementation)) implementation <- "cpp"

    if(implementation == "cpp"){

      result <- csmsb_allocate_cpp(data = data,
                                   arm = arm,
                                   newData = newData,
                                   settings = settings,
                                   randSequence = randSequence
      )

    } else if (implementation == "R"){

      result <- csmsb_allocate_R(data = data,
                                   arm = arm,
                                   newData = newData,
                                   settings = settings,
                                   randSequence = randSequence
      )

    }

  } else {
    stop("Under construction")
  }

  result

}


