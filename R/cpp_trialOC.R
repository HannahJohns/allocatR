#' Operating Characteristics
#'
#' @description Estimate the operating characteristics of a trial using the MSB, CS-MSB or CSSize-MSB algorithms
#'
#' @usage
#' trialOC(data, outcomes=NULL, primaryAnalysis = NULL,
#'          n=NULL, replay=FALSE,
#'          alpha=0.05,
#'          method='CSMSB',
#'          voteMethod = NULL,
#'          biasMethod = "static",
#'          weights=NULL,
#'          biasParameters,
#'          groupImbalanceThreshold=Inf,
#'          imbalancePriority="abstain",
#'          nThreads=NULL,
#'          nSims,
#'          trng = NULL,
#'          jump_trng = 0
#' )
#'
#' @param data A data frame of covariates to sample from
#' @param outcomes An object containing counterfactual outcome information for each record. See details.
#' @param primaryAnalysis A user generated function which performs analysis on the simulated data
#' @param n If replay is false, the sample size for the trial
#' @param replay If the trial should recruit exactly who is specified in data (implies n=nrow(data))
#' @param alpha A number indicating the level of significance for when MSB/CSMSB intervenes.
#' @param method A character string indicating the type of method to use.  Defaults to CSMSB.
#' @param voteMethod A character string indicating how votes should be combined. Defaults to "majority" for MSB and "pooled" for CSMSB.
#' @param biasMethod A character string indicating how a bias coin probability should be generated. Defaults to "static"
#' @param weights A vector of length \code{ncol(data)} giving the weights for each column in order. Defaults to even weights if not supplied.
#' @param biasParameters A numeric vector providing parameters to the indicated bias method
#' @param groupImbalanceThreshold Group allocation imbalance threshold required to trigger CSSize-MSB allocation rules. By default will never trigger.
#' @param imbalancePriority Priority rule for CSSize-MSB. Defaults to "abstain"
#' @param nThreads Number of parallel threads to use in simulation. Defaults to maximum available if the parallel package is installed, otherwise 1
#' @param nSims Number of trials to simulate to estimate operating characteristics
#' @param trng an optional parallel random number generator from the rTRNG package. If an integer is supplied instead, a yarn3 engine will be generated using this integer as the seed.
#' @param jump_trng a nonnegative integer specifying how many steps the trng engine should be jumped forward before it is used.
#'
#'
#' @returns A data frame of nSims rows containing the following columns:
#'
#' \itemize{
#' \item dataHash - An md5 hash of the input data
#' \item method - The specified method
#' \item simID - The simulation number
#' \item n - The specified sample size. NA if not supplied.
#' \item alpha - The specified alpha
#' \item voteMethod - The specified vote method
#' \item biasMethod - The specified bias method
#' \item biasParamters - The specified bias parameters
#' \item replay - The specified value for replay
#' \item expected_bias - The mean deviation from a fair coin
#' \item prop_random - The proportion of randomisations using a fair coin
#' \item prop_random_covariate  - The proportion of randomisations where a fair coin would be used for covariate balance
#' \item prop_random_group  - The proportion of randomisations where a fair coin would be used for group imbalance balance
#' \item std_imbalance_... - Test statistics for imbalance at end of trial using t tests and chi-square tests
#' \item std_p_... - Test statistics for imbalance at end of trial using p values (t tests and chi-square tests)
#' \item genodds_imbalance_... Test statistics for imbalance at end of trial using genodds
#' \item genodds_p_.... Test statistics for imbalance at end of trial using p values (genodds)
#' \item vote_vote0_... Proportion of times each variable voted to assign to control
#' \item vote_vote1_...  Proportion of times each variable voted to assign to treatment
#' \item vote_abstain_... Proportion of times each variable voted to abstain
#' }
#'
#'
trialOC <- function(data,
                    outcomes=NULL,
                    primaryAnalysis = NULL,
                    n=NULL,
                    replay=FALSE,
                    alpha=0.05,
                    method='CSMSB',
                    voteMethod = NULL,
                    biasMethod = "static",
                    weights=NULL,
                    biasParameters,
                    groupImbalanceThreshold=Inf,
                    imbalancePriority="abstain",
                    nThreads=NULL,
                    nSims,
                    trng = NULL,
                    jump_trng = 0
                    ){

  if(method == "CSMSB"){
    if(is.null(voteMethod)) voteMethod <- "pooled"
  } else if(method=="MSB"){
    if(is.null(voteMethod)) voteMethod <- "majority"
    if(voteMethod != "majority") stop("MSB algorithm only supports majority voting")
  }

  if(voteMethod == "pooled"){
    if(is.null(weights)){
      weights <- rep(1,ncol(data))
    }
  }

  if(!is.null(weights)){
    if(voteMethod != "pooled") top("weights are only supported if voteMethod is pooled")
    if(length(weights) != ncol(data)) stop("length of weights should match columns in data")
  }


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

  # Run a single simulation and compress the results to a single line
  batchSimFunction <- function(i,data,outcomes,primaryAnalysis,n,
                               replay,alpha,method,voteMethod, weights,
                               biasMethod,biasParameters,
                               groupImbalanceThreshold,imbalancePriority,
                               trng, jump_trng
                              ){

    # If trng is used, need to calculate how many random numbers are used in each simulation.


    # If replay is false, we use n when resampling to get the data for the simulation run.
    # We need either n or nrow(data) random uniform numbers per run
    # Once outcome generation is enabled, this will need to go in as well.

    if(replay){
      num_sims_per_run <- nrow(data)
    } else {
      if(is.null(n)) stop("n must be supplied if replay = FALSE")
      num_sims_per_run <- 2*n
    }

    # Get a hash of the data. Used to ensure comparisons are made on the same data
    dataHash <- digest::digest(data,algo = "md5")

    thisSim <- csmsb:::trialSim(
                        data = data,
                        outcomes = outcomes,
                        primaryAnalysis = primaryAnalysis,
                        n = n,
                        replay = replay,
                        alpha = alpha,
                        method = method,
                        voteMethod = voteMethod,
                        biasMethod = biasMethod,
                        weights=weights,
                        biasParameters = biasParameters,
                        groupImbalanceThreshold=groupImbalanceThreshold,
                        imbalancePriority=imbalancePriority,
                        trng = trng,
                        jump_trng = jump_trng+(i-1)*num_sims_per_run
                        )

    # Flatten the results to a single row
    expectedAbsoluteBias <-  thisSim$expected_bias
    prop_random  <- thisSim$prop_random
    prop_random_covariate  <- thisSim$prop_random_covariate
    prop_random_group  <- thisSim$prop_random_group

    imbalance_std <- thisSim$imbalance_std
    names(imbalance_std) <- paste("std",names(imbalance_std),sep="_")

    imbalance_genodds <- thisSim$imbalance_genodds
    names(imbalance_genodds) <- paste("genodds",names(imbalance_genodds),sep="_")

    voteDist <- c(thisSim$vote_dist)
    names(voteDist) <- c(outer(rownames(thisSim$vote_dist),colnames(thisSim$vote_dist),paste,sep="_"))
    names(voteDist) <- paste("vote",names(voteDist),sep="_")

    if(is.null(n)) n <- NA
    sim_results <- data.frame(
                     simID=i,
                     n=n,
                     alpha=alpha,
                     voteMethod = voteMethod,
                     biasMethod = biasMethod,
                     biasParameters=paste(biasParameters,collapse=" "),
                     replay=replay,
                     expected_bias=expectedAbsoluteBias,
                     prop_random=prop_random,
                     prop_random_covariate=prop_random_covariate,
                     prop_random_group=prop_random_group,
                     n0=sum(thisSim$arm=="control"),
                     n1=sum(thisSim$arm=="treatment")
                     )

    sim_results$allocation_difference = abs(sum(thisSim$arm=="control") - sum(thisSim$arm=="treatment"))

    sim_results <- cbind(sim_results,t(imbalance_std),t(imbalance_genodds),t(voteDist))
    sim_results <- cbind(dataHash=dataHash,method=method,sim_results)

    if(!is.null(primaryAnalysis))
    {
      sim_results <- cbind(sim_results,
                           as.data.frame(t(as.matrix(thisSim$analysisResults)))
      )
    }

    return(sim_results)
  }

  # Either wrap snow functions or create dummies
  if(nThreads>1)
  {
    cl <- parallel::makeCluster(nThreads)
    batchRun <- parallel::parLapply(cl, 1:nSims, batchSimFunction,
                                    data=data,
                                    outcomes=outcomes,
                                    primaryAnalysis=primaryAnalysis,
                                    n=n,
                                    replay=replay,
                                    alpha=alpha,
                                    method=method,
                                    voteMethod=voteMethod,
                                    biasMethod=biasMethod,
                                    weights=weights,
                                    groupImbalanceThreshold=groupImbalanceThreshold,
                                    imbalancePriority=imbalancePriority,
                                    biasParameters=biasParameters,
                                    trng=trng,
                                    jump_trng=jump_trng
    )
    parallel::stopCluster(cl)

  } else {

    batchRun <- lapply(1:nSims, batchSimFunction,
                       data=data,
                       outcomes=outcomes,
                       primaryAnalysis=primaryAnalysis,
                       n=n,
                       replay=replay,
                       alpha=alpha,
                       method=method,
                       voteMethod=voteMethod,
                       biasMethod=biasMethod,
                       weights=weights,
                       groupImbalanceThreshold=groupImbalanceThreshold,
                       imbalancePriority=imbalancePriority,
                       biasParameters=biasParameters,
                       trng=trng,
                       jump_trng=jump_trng
    )
  }

  batchRun <- do.call("rbind",batchRun)

  return(batchRun)

}
