trialSim <- function(  method = c("simple","csmsb","msb"),
                       data,
                       replay = FALSE,
                       settings = NULL,
                       randSequence = NULL,
                       implementation = NULL,
                       simple=NULL,
                       nThreads = NULL
){

  # Generic input sanitation
  data <- as.data.frame(data)

  #' trialSim() runs recursively if multiple randSequences are provided

  if("list" %in% class(randSequence$sequence)){

    if(is.null(simple)){
      simple <- TRUE
    }

    ## Recursive Call ########################
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
                                 method,
                                 data,
                                 replay,
                                 settings,
                                 randSequence,
                                 implementation,
                                 simple
    ){

      # All sequences should be the same size

      # We need to know total sample size because we
      # mihgt randomly resample rows in each sim.
      # If we do, we need toto jump the seed state ahead
      # the correct ammount to avoid double-dipping

      if("Runif" %in% class(randSequence$sequence[[1]])) {
        n <- length(randSequence$sequence[[1]])
      } else if("Permutedblock" %in% class(randSequence$sequence[[1]])) {
        n <- nrow(randSequence$sequence[[1]])
      } else {
        stop("Unrecognised sequence type")
      }

      if(replay){
        num_to_jump_per_run <- 0 # We're not resampling rows, so there's no extra randomness
      } else {
        num_to_jump_per_run <- n # We need n random numbers per run
      }

      # Use a modified randSeq containing just the i-th sequence plus appropriate
      # RNG engine

      this_sequence = randSequence$sequence[[i]]
      this_trng = randSequence$trng$copy()
      this_trng$jump((i-1) * num_to_jump_per_run)

      this_randSeq <- list(
       sequence = this_sequence,
       trng = this_trng
      )

      thisSim <- trialSim(method = method,
                          data=data,
                          replay = replay,
                          settings = settings,
                          randSequence = this_randSeq,
                          implementation = implementation,
                          simple=simple
      )

      return(thisSim)
    }




    # Either wrap snow functions or create dummies
    if(nThreads>1)
    {
      cl <- parallel::makeCluster(nThreads)

      batchRun <- parallel::parLapply(cl, 1:nSims, batchSimFunction,
                                      method=method,
                                      data=data,
                                      replay=replay,
                                      settings=settings,
                                      randSequence=randSequence,
                                      implementation=implementation,
                                      simple=simple
      )
      parallel::stopCluster(cl)

    } else {

      batchRun <- lapply(1:nSims, batchSimFunction,
                         method=method,
                         data=data,
                         replay=replay,
                         settings=settings,
                         randSequence=randSequence,
                         implementation=implementation,
                         simple=simple
      )

    }

    if(simple){
      batchRun <- do.call("rbind",batchRun)
    }

    return(batchRun)

  } else {

    if(!("Alloc_seq" %in% class(randSequence$sequence))) stop("randSequence should contain an allocation sequence")

    ## Base Call ########################

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

    #' NOTE: we create a copy of the engine, rather than passing by reference.
    #' This means that if a user runs trialSim() multiple times using a saved
    #' random sequence, they will get the same resampled data each time as well.
    #'
    #' This is (probably) the correct behaviour to use here because it's internally
    #' consistent and ensures reproducibility of results.
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

      if(is.null(implementation)) implementation <- "cpp"

      if(implementation == "cpp"){

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


      } else if (implementation == "R"){
        stop("Under construction")

      } else stop("unrecognised implementation")

    } else {
      stop("Unrecognised method")
    }

    return(result)

  }


}
