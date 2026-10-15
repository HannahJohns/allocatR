#' Randomise a participant
#'
#' @description Randomise a participant.
#'
#' @usage allocate( randomisationScheme,
#'                  randSequence,
#'                  arm,
#'                  data = NULL,
#'                  newData = NULL,
#'                  offset=0,
#'                  implementation = c("cpp", "R")
#'                  )
#'
#' @param randomisationScheme an object of type \link{RandomisationScheme}
#' @param arm An integer vector containing the current treatment assignments (0=control, 1=intervention)
#' @param data A data frame with \code{length(arm)} rows containing relevant data (e.g. covariates) for currently assigned participants
#' @param newData A one-row data frame with the same columns as \code{data} containning relevant data (e.g. covariates) for the new participant to be assigned
#' @param offset A positive integer to increase the allocation sequence step by. Used to prevent re-use of the allocation sequence.
#' @param randSequence An object of type \link{Alloc_seq}
#' @param implementation A string indicating which implementation of a given method should be used
#'
#' @returns A list with class \code{Allocation} containing the following:
#' \describe{
#' \item assignment
#' \item details a list containing internal steps taken to reach this assignment result
#' }
#'
#' @details
#'
#' implementation should be left to the default (cpp) unless specifically required due to
#' e.g. system architecture/OS constraints.
#'
#' If relevant, variable types are read from the class of each column in \code{data}.
#' Factors are assumed categorical, ordered factors are assumed ordinal and numeric
#' are assumed continuous.
#'
#' See details of functions that produce \link[=RandomisationScheme]{randomisationScheme} for specific details
#' of each method.
#'
#' @examples
#'
#' data <- synthetic_covariate_data[1:30,]
#' newData <- synthetic_covariate_data[31,,drop=F]
#' arm <- rep(0:1, length.out = 30)
#'
#' randScheme <- csmsb(alpha = 0.5, biasParameters=0.7)
#' randSeq <- get_randNumbers(50)
#'
#' allocate(randomisationScheme=randScheme,  randSequence=randSeq, arm=arm, data=data, newData=newData)
#'
allocate <- function(  randomisationScheme,
                       randSequence,
                       arm,
                       data = NULL,
                       newData = NULL,
                       offset = 0,
                       implementation = NULL
){

  method <- randomisationScheme$method
  settings <-  randomisationScheme$settings

  # Generic input sanitation
  data <- as.data.frame(data)
  newData <- as.data.frame(newData)

  if(method == "csmsb"){

    if(is.null(implementation)) implementation <- "cpp"

    if(implementation == "cpp"){

      result <- allocatR::csmsb_allocate_cpp(data = data,
                                   arm = arm,
                                   newData = newData,
                                   settings = settings,
                                   randSequence = randSequence,
                                   offset = offset
      )

    } else if (implementation == "R"){

      result <- allocatR::csmsb_allocate_R(data = data,
                                   arm = arm,
                                   newData = newData,
                                   settings = settings,
                                   randSequence = randSequence,
                                   offset = offset
      )

    }

  } else if(method == "simple"){

    if(!is.null(implementation)) warning("implementation option ignored for simple randomisation scheme")

    result <- allocatR::simple_allocate(data = data,
                               arm = arm,
                               newData = newData,
                               settings = settings,
                               randSequence = randSequence,
                               offset = offset
    )


  } else if(method == "msb") { 
  
    if(is.null(implementation)) implementation <- "cpp"

    if(implementation == "cpp"){

      result <- allocatR::msb_allocate_cpp(data = data,
                                           arm = arm,
                                           newData = newData,
                                           settings = settings,
                                           randSequence = randSequence,
                                           offset = offset
      )

    } else if (implementation == "R"){

      result <- allocatR::msb_allocate_R(data = data,
                                         arm = arm,
                                         newData = newData,
                                         settings = settings,
                                         randSequence = randSequence,
                                         offset = offset
      )

    }

  } else {
    stop("unrecognised randomisation scheme")
  }

  result

}