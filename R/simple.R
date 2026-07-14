#' Simple Randomisation
#'
#' @description
#' Creates a randomisation scheme using simple randomisation
#'
#' @usage
#' simple()
#'
#' @details
#' Performs either simple randomisation or permuted block randomisation with no
#' adaptive components.
#'
#' @returns an object of class RandomisationScheme
#'
simple <- function(){
  out <- list(
    method = "simple",
    settings = NULL
  )

  class(out) <- c("RandomisationScheme",class(out))
  out
}


simple_allocate <- function(data,
                            arm,
                            newData,
                            settings,
                            randSequence
){
  # Extract out and process inputs
  rand_number <- length(arm)+1


  # Check that only one sequence provided
  if("list" %in% class(randSequence$sequence)){
    if("Alloc_seq" %in% class(randSequence$sequence[[1]])){
      stop("Only one allocation sequence should be provided, not a list of them.")
    }
  }

  # Get the random number for this allocation

  if("Runif" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence
    roll <- randSequence[rand_number]

  } else if ("Permutedblock" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence
    roll <- randSequence[rand_number,1]
  }


  decision <- ifelse(roll<x$probability,0,1)

  return(list(
    decision = decision,
    prob = 0.5,
    roll = roll,
    details = NULL
  ))

}



simple_trialSim <- function(data,
                            settings,
                            randSequence,
                            replay
){

  # Generate random numbers for randomisation and feed this list to C++
  if("Runif" %in% class(randSequence$sequence)){

    roll <- randSequence$sequence

  } else if ("Permutedblock" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence

    roll <- randSequence[,1]

  }

  arm <- ifelse(roll<0.5, 0, 1)

  # Format results


  arm <- factor(arm,levels=0:1,labels=c("control","treatment"))

  probabilityLog <- probabilityLog_covariate <- probabilityLog_group <- rep(0.5, length(arm))

  out <- list(probabilityLog=probabilityLog,
              probabilityLog_covariate=probabilityLog_covariate,
              probabilityLog_group=probabilityLog_group,
              rollLog=roll,
              arm = arm)

  out$imbalance_std <- getImbalance(data,arm,standard=T)
  out$imbalance_genodds <- getImbalance(data,arm,standard = F)
  out$prop_random <- getProportionRandom(probabilityLog)

  out$prop_random_covariate <- getProportionRandom(probabilityLog_covariate)
  out$prop_random_group <- getProportionRandom(probabilityLog_group)

  out$expected_bias <- getExpectedBias(probabilityLog)

  return(out)

}
