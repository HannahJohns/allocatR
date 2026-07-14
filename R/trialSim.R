#' Simulate trials
#'
#' @description
#' Randomise a participant according to the MSB, CS-MSB or CSSize-MSB algorithms.
#'
#' @usage trialSim(randomisationScheme,
#'                 randSequence = NULL,
#'                 data = NULL,
#'                 replay = FALSE,
#'                 simple=NULL,
#'                 nThreads=NULL
#'                 )
#'
#' @param randomisationScheme an object of type \link{RandomisationScheme}
#' @param randSequence An object of type \link{Alloc_seq}
#' @param data A data frame containing relevant data (e.g. covariates) for participants in simulated trials
#' @param replay A logical value. If \code{TRUE} then trial participants are randomised in sequence from the first to last row of \code{data}.
#'               If \code{FALSE}, then participants are randomly sampled with replacement from \code{data}
#' @param simple A logical value. If \code{TRUE} then results of each simulated trial are condensed to a single row.
#' @param nThreads Number of parallel threads to use in simulation.
#'
#'
#' @details
#'
#' Size of trial and number of simulated trials are both specified by \link[=Alloc_seq]{randSequence} object.
#'
#' If \code{replay} is \code{TRUE} then the trial length must equal \code{nrow(data)}.
#'
#' If relevant, variable types are read from the class of each column in \code{data}.
#' Factors are assumed categorical, ordered factors are assumed ordinal and numeric
#' are assumed continuous.
#'
#' See details of functions that produce \link[=RandomisationScheme]{randomisationScheme} for specific details
#' of each method.

#' @returns
#'
#' If \code{simple} is \code{FALSE}, returns a list of results for the simulated trial. This includes:
#'
#' The following participant-level details:
#' \describe{
#' \item probabilityLog - the probability of assignment to control for each participant
#' \item probabilityLog_covariate - the probability of assignment to control for each participant, just based on covariate information
#' \item probabilityLog_group - the probability of assignment to control for each participant, just based on group allocation information alone
#' \item rollLog - the random number used to determine treatment assignment
#' \item arm - the resulting treatment allocation
#' }
#'
#' And the following participant-level details:
#' \describe{
#' \item imbalance_std - covariate balance measured using standard tests (test statistic and p-value from the t-test or chi-square test)
#' \item imbalance_genodds - covariate balance measured using generalised odds ratios (WMW OR and p-value. Categorical variables are one-hot coded for this evaluation)
#' \item prop_random - proportion of interventions where probability of assignment was 0.5
#' \item prop_random_covariate - proportion of interventions where probability of assignment (based on covariate information alone) was 0.5
#' \item prop_random_group - proportion of interventions where probability of assignment (based on group allocaiton information alone) was 0.5
#' \item expected_bias - the average deviation from a probability of 0.5 across all group allocations.
#' }
#'
#' Plus any specific outputs relevant to randomisationScheme.
#'
#' If \code{simple} is \code{TRUE}, returns a one-row list of summary results containing just the trial-level details above.
#'
#' If more than one trial was run, The above is returned either as a list (if \code{simple} is \code{FALSE}) or a data frame (if \code{simple} is \code{TRUE})
#'
#' @examples
#' data <- synthetic_covariate_data
#'
#' randScheme <- csmsb(alpha = 0.5, biasParameters=0.7)
#' randSeq <- get_randNumbers(100, trng=1)
#' trialSim(randomisationScheme=randScheme, randSequence=randSeq,  data=data, replay=FALSE)
#'
#' randSeq2 <- get_randNumbers(100, n.trials=10, trng=2)
#' trialSim(randomisationScheme=randScheme, randSequence=randSeq2,  data=data, replay=FALSE)
trialSim <- function(randomisationScheme,
                     randSequence,
                     data = NULL,
                     replay = FALSE,
                     simple=NULL,
                     nThreads = NULL
){


  # This errors when called recursively
  if(!("RandomisationScheme" %in% class(randomisationScheme))) stop("randomisationScheme should be a randomisation scheme")
  if(!("Alloc_seq" %in% class(randSequence$sequence))){
    if("list" %in% class(randSequence$sequence)){
      if(!("Alloc_seq" %in% class(randSequence$sequence[[1]]))){
        stop("randSequence not properly formatted")
      }
    } else{
      stop("randSequence not properly formatted (should contain list or Alloc_seq object)")
    }
  }

  # stop("We're losing class information from randSequence for some reason?")
  method <- randomisationScheme$method
  settings <-  randomisationScheme$settings

  # Generic input sanitation
  data <- as.data.frame(data)

  # trialSim() runs recursively if multiple randSequences are provided
  if("list" %in% class(randSequence$sequence)){

    ## Recursive Call ########################


    if(is.null(simple)){
      simple <- TRUE
    }

    nSims <- length(randSequence$sequence)
    if(is.null(nThreads)){
      if("parallel" %in% installed.packages()[,1])
      {
        nThreads <- parallel::detectCores()
      } else {
        nThreads <- 1
      }
    }

    if(nThreads>1 & !("parallel" %in% installed.packages()[,1])){
      warning("multiThreadding requires the parallel package. Please run install.packages(\"parallel\")")
      nThreads <- 1
    }

    # Run a single simulation and compress the results to a single line.
    # This is just a fancy wrapper for trialSim, which gets called recursively.

    batchSimFunction <- function(i,
                                 randomisationScheme,
                                 data,
                                 replay,
                                 randSequence_details,
                                 simple
    ){

      this_sequence = randSequence_details$sequence[[i]]

      this_trng <- switch(randSequence_details$trng_type,
                          Rcpp_yarn2=rTRNG::yarn2$new(randSequence_details$trng_string),
                          Rcpp_yarn3=rTRNG::yarn3$new(randSequence_details$trng_string),
                          Rcpp_yarn4=rTRNG::yarn4$new(randSequence_details$trng_string),
                          Rcpp_yarn5=rTRNG::yarn5$new(randSequence_details$trng_string),
                          Rcpp_yarn3s=rTRNG::yarn3s$new(randSequence_details$trng_string),
                          Rcpp_yarn5s=rTRNG::yarn5s$new(randSequence_details$trng_string),
                          Rcpp_lcg64=rTRNG::lcg64$new(randSequence_details$trng_string),
                          Rcpp_lcg64_shift=rTRNG::lcg64_shift$new(randSequence_details$trng_string),
                          Rcpp_mrg2=rTRNG::mrg2$new(randSequence_details$trng_string),
                          Rcpp_mrg3=rTRNG::mrg3$new(randSequence_details$trng_string),
                          Rcpp_mrg4=rTRNG::mrg4$new(randSequence_details$trng_string),
                          Rcpp_mrg5=rTRNG::mrg5$new(randSequence_details$trng_string),
                          Rcpp_mrg3s=rTRNG::mrg3s$new(randSequence_details$trng_string),
                          Rcpp_mrg5s=rTRNG::mrg5s$new(randSequence_details$trng_string),
                          stop("UNRECOGNISED TRNG ENGINE TYPE")
      )

      # All sequences should be the same size

      # We need to know total sample size because we
      # mihgt randomly resample rows in each sim.
      # If we do, we need toto jump the seed state ahead
      # the correct ammount to avoid double-dipping


      if("Runif" %in% class(this_sequence)) {
        n <- length(this_sequence)
      } else if("Permutedblock" %in% class(this_sequence)) {
        n <- nrow(this_sequence)
      } else {
        stop("Unrecognised sequence type")
      }

      if(replay){
        num_to_jump_per_run <- 0 # We're not resampling rows, so there's no extra randomness
      } else {
        num_to_jump_per_run <- n # We need n randthis_sequenceom numbers per run
      }

      # Use a modified randSeq containing just the i-th sequence plus appropriate
      # RNG engine
      this_trng$jump((i-1) * num_to_jump_per_run)

      this_randSeq <- list(
        sequence = this_sequence,
        trng = this_trng
      )

      thisSim <- trialSim(randomisationScheme=randomisationScheme,
                          data=data,
                          replay = replay,
                          randSequence = this_randSeq,
                          simple=simple
      )

      return(thisSim)
    }

    # Feeding through randSequence directly causes problems. Instead, feed
    # through enough detail to reconstruct it inside the loop
    randSequence_details = list(
      sequence = randSequence$sequence,
      trng_type = class(randSequence$trng),
      trng_string = randSequence$trng$toString(),
      class = class(randSeq2$sequence[[1]])
    )

    # Either wrap snow functions or create dummies
    if(nThreads>1)
    {
      cl <- parallel::makeCluster(nThreads)

      batchRun <- parallel::parLapply(cl, 1:nSims, batchSimFunction,
                                      randomisationScheme=randomisationScheme,
                                      data=data,
                                      replay=replay,
                                      randSequence_details=randSequence_details,
                                      simple=simple
      )
      parallel::stopCluster(cl)

    } else {

      batchRun <- lapply(1:nSims, batchSimFunction,
                         randomisationScheme=randomisationScheme,
                         data=data,
                         replay=replay,
                         randSequence_details=randSequence_details,
                         simple=simple
      )

    }

    if(simple){
      batchRun <- do.call("rbind",batchRun)
    }

    return(batchRun)

  } else {

    ## Base Case ########################

    if(!("Alloc_seq" %in% class(randSequence$sequence))) stop("randSequence should contain an allocation sequence")

    if(is.null(simple)) simple <- FALSE

    # Get the sample size
    if("Runif" %in% class(randSequence$sequence)) {
      n <- length(randSequence$sequence)
    } else if("Permutedblock" %in% class(randSequence$sequence)) {
      n <- nrow(randSequence$sequence)
    } else {
      stop("Unrecognised sequence type")
    }


    # Resample data if this was specified
    # NOTE: we create a copy of the engine, rather than passing by reference.
    # This means that if a user runs trialSim() multiple times using a saved
    # random sequence, they will get the same resampled data each time as well.
    #
    # This is (probably) the correct behaviour to use here because it's internally
    # consistent and ensures reproducibility of results.
    trng <- randSequence$trng$copy()
    rownames(data) <- NULL
    if(replay){
      if(n != nrow(data)) stop("Allocation sequence should be of length nrow(data) if replay is TRUE")
    } else {
      # If we're not replaying the dataset, we need to sample rows to feed through
      if(is.null(trng)){
        rows <- ceiling(runif(n,min=0,max=nrow(data)))
      } else {
        rows <- ceiling(rTRNG::runif_trng(n,min=0,max=nrow(data),engine = trng))
      }

      data <- data[rows,]

    }

    ### Method specific calls ########################

    if(method == "csmsb"){

      result <- csmsb_trialSim_cpp(data = data,
                                   settings = settings,
                                   randSequence = randSequence,
                                   replay=replay
      )

      if(simple){

        # Flatten the results to a single row


        expectedAbsoluteBias <-  result$expected_bias
        prop_random  <- result$prop_random
        prop_random_covariate  <- result$prop_random_covariate
        prop_random_group  <- result$prop_random_group

        imbalance_std <- result$imbalance_std
        names(imbalance_std) <- paste("std",names(imbalance_std),sep="_")

        imbalance_genodds <- result$imbalance_genodds
        names(imbalance_genodds) <- paste("genodds",names(imbalance_genodds),sep="_")

        voteDist <- c(result$vote_dist)
        names(voteDist) <- c(outer(rownames(result$vote_dist),colnames(result$vote_dist),paste,sep="_"))
        names(voteDist) <- paste("vote",names(voteDist),sep="_")


        for(i in 1:length(settings)){
          if(is.null(settings[[i]])) settings[[i]] <- NA
          if(length(settings[[i]]) > 1) settings[[i]] <- paste(settings[[i]], collapse = ", ")
        }


        sim_results <- cbind(
          do.call("data.frame",settings),
          data.frame(
            n=n,
            expected_bias=expectedAbsoluteBias,
            prop_random=prop_random,
            prop_random_covariate=prop_random_covariate,
            prop_random_group=prop_random_group,
            n0=sum(result$arm=="control"),
            n1=sum(result$arm=="treatment")
          )

        )

        sim_results$allocation_difference = abs(sum(result$arm=="control") - sum(result$arm=="treatment"))

        sim_results <- cbind(sim_results,t(imbalance_std),t(imbalance_genodds),t(voteDist))
        sim_results <- cbind(method=method,sim_results)

        result <- sim_results

      }

    } else if(method == "simple"){

      result <- allocatR::simple_trialSim(data = data,
                                          settings = settings,
                                          randSequence = randSequence,
                                          replay=replay
      )

      if(simple){

        # Flatten the results to a single row

        expectedAbsoluteBias <-  result$expected_bias
        prop_random  <- result$prop_random
        prop_random_covariate  <- result$prop_random_covariate
        prop_random_group  <- result$prop_random_group

        imbalance_std <- result$imbalance_std
        names(imbalance_std) <- paste("std",names(imbalance_std),sep="_")

        imbalance_genodds <- result$imbalance_genodds
        names(imbalance_genodds) <- paste("genodds",names(imbalance_genodds),sep="_")


        sim_results <- cbind(
          data.frame(
            n=n,
            expected_bias=expectedAbsoluteBias,
            prop_random=prop_random,
            prop_random_covariate=prop_random_covariate,
            prop_random_group=prop_random_group,
            n0=sum(result$arm=="control"),
            n1=sum(result$arm=="treatment")
          )

        )

        sim_results$allocation_difference = abs(sum(result$arm=="control") - sum(result$arm=="treatment"))

        sim_results <- cbind(sim_results,t(imbalance_std),t(imbalance_genodds))
        sim_results <- cbind(method=method,sim_results)

        result <- sim_results

      }


    } else {
      stop("Unrecognised method")
    }

    return(result)


  }

  return(NULL)
}
