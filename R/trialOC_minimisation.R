trialOC_minimisation <- function(data,stratum=NULL,
                    n=NULL,
                    replay=FALSE,
                    biasStrength,
                    balanceGroups=FALSE,
                    nThreads=NULL,
                    rollBlocks=NULL,
                    interrims=NULL,
                    nSims){

  if(is.null(nThreads))
  {
    if("parallel" %in% installed.packages()[,1])
    {
      nThreads <- parallel::detectCores()
    } else {
      nThreads <- 1
    }
  }

  if(nThreads>1 & !("parallel" %in% installed.packages()[,1]))
  {
    warning("multiThreadding requires the parallel package. Please run install.packages(\"parallel\")")
    nThreads <- 1
  }

  # Run a single simulation and compress the results to a single line
  batchSimFunction <- function(i,data,stratum,outcomes,n,
                               replay,biasStrength,balanceGroups,interrims
  ){

    # Get a hash of the data. Used to ensure comparisons are made on the same data
    dataHash <- digest::digest(data,algo = "md5")

    thisSim <- csmsb:::trialSim_minimisation(data=data,stratum=stratum,
                                n=n,
                                replay=replay,
                                biasStrength=biasStrength,
                                balanceGroups = balanceGroups,
                                rollBlocks=rollBlocks,
                                interrims=interrims)


    # Get n for the results if we were using sample size
    if(replay) n <- nrow(data)

    # Flatten the results to a single row
    expectedAbsoluteBias <-  thisSim$expected_bias
    prop_random  <- thisSim$prop_random

    imbalance_std <- thisSim$imbalance_std
    colnames(imbalance_std) <- paste("std",colnames(imbalance_std),sep="_")

    imbalance_genodds <- thisSim$imbalance_genodds
    colnames(imbalance_genodds) <- paste("genodds",colnames(imbalance_genodds),sep="_")

    # voteDist <- c(thisSim$vote_dist)
    # names(voteDist) <- c(outer(rownames(thisSim$vote_dist),colnames(thisSim$vote_dist),paste,sep="_"))
    # names(voteDist) <- paste("vote",names(voteDist),sep="_")

    allocationBalance <- thisSim$allocationBalance
    colnames(allocationBalance) <- paste("n",colnames(allocationBalance),sep="_")

    sim_results <- cbind(simID=i,
                     n=n,
                     rollUniformCheck=thisSim$rollUniformCheck,
                     biasStrength=biasStrength,
                     replay=replay,
                     look=thisSim$look,
                     stratum_balance=thisSim$stratum_balance,
                     expected_bias=expectedAbsoluteBias,
                     prop_random=prop_random,
                     allocationBalance,
                     imbalance_std,
                     imbalance_genodds
                     )
    sim_results <- as.data.frame(sim_results)
    sim_results <- cbind(dataHash=dataHash,sim_results)

    return(sim_results)
  }

  # Either wrap snow functions or create dummies
  if(nThreads>1)
  {
    cl <- parallel::makeCluster(nThreads)
    batchRun <- parallel::parLapply(cl, 1:nSims, batchSimFunction,
                                    data=data,
                                    stratum=stratum,
                                    n=n,
                                    replay=replay,
                                    biasStrength=biasStrength,
                                    balanceGroups=balanceGroups,
                                    interrims=interrims
    )
    parallel::stopCluster(cl)

  } else {

    batchRun <- lapply(1:nSims, batchSimFunction,
                       data=data,
                       stratum=stratum,
                       n=n,
                       replay=replay,
                       biasStrength=biasStrength,
                       balanceGroups=balanceGroups,
                       interrims=interrims
    )
  }

  batchRun <- do.call("rbind",batchRun)

  return(batchRun)

}
