#' Common Scale Minimal Sufficient Balance
#'
#' @description
#' Creates a randomisation scheme for the CSMSB or CS-Size MSB algorithms
#'
#' @usage
#' csmsb(alpha,
#'       biasParameters,
#'       voteMethod="pooled",
#'       biasMethod="static",
#'       weights = NULL,
#'       groupImbalanceThreshold=Inf,
#'       imbalancePriority = "abstain"
#' )
#'
#' @param alpha A number indicating the level of significance for when CSMSB/CS-Size intervenes based on covariate imbalance.
#' @param biasParameters A numeric vector providing parameters to the indicated bias method. See details.
#' @param voteMethod A character string indicating how covariate votes should be combined.
#' @param biasMethod A character string indicating how a bias coin probabilities should be generated.
#' @param weights A named vector giving the weights for each covariate. Defaults to even weights if no value is supplied.
#' @param groupImbalanceThreshold Group allocation imbalance threshold required to trigger CSSize-MSB allocation rules. By default this will never trigger (i.e. CSMSB is used, not CSSize-MSB).
#' @param imbalancePriority Priority rule for CSSize-MSB. Defaults to "abstain".
#'
#' @details
#'
#' \code{biasParameters} should be either a single number or a vector of numbers
#' depending on \code{biasMethod}
#'
#' \enumerate{
#' \item static - a single number between 0.5 and 1 giving the strength of the biased coin
#' \item linear - a single number between 1 and infinity that scales the imbalance measure to calculate a biased coin strength according to the formula 0.5+imbalance*biasParameters
#' }
#'
#' When csmsb() is used with the \link[=trialSim]{trialSim} function and with majority or weighted
#' voting, there are a few quirks to its implementation.
#' \enumerate{
#'     \item When checking each variable, we assume that a vote is always valid and log the vote as if it is.
#'     \item We remove invalid votes based on \code{alpha}, immediately before we tally votes.
#' }
#' This is mainly done so we can log small imbalances for understanding trial behaviour over time.
#'
#' When csmsb() is used with the \link[=allocate]{allocate} function, the details list in its output
#' will contain the following:
#'   \itemize{
#'            \item pVals - p-values indicating imbalance for each covariate
#'            \item vote - The votes placed by each covariate
#'            \item imbalance - The current imbalance on each covariate
#'            \item seImbalance - The standard error of the imbalance on each covariate
#'            \item prospectiveImbalance0 - The imalance on each covariate if the new participant is assigned to control,
#'            \item prospectiveImbalance1 - The imalance on each covariate if the new participant is assigned to intervention,
#'            \item BD - The pooled measure of imbalance,
#'            \item seBD - The standard error of the pooled measure of imbalance,
#'            \item BD_pVal - The p-value of the pooled measure of imbalance,
#'            \item covBalanceProb=covBalanceProb - the probability of assignent to control if only covariate balance is considered,
#'            \item groupImbalanceProb=groupImbalanceProb  - the probability of assignent to control if only group allocation balance is considered
#'     }
#'
#'
#' When csmsb() is used with the \link[=trialSim]{trialSim} function, the details list in its output
#' will additionally contain the following:
#'   \itemize{
#'            \item voteLog - the votes placed by each covariate during each allocation
#'            \item pValLog - covariate-level p-value for current imbalance
#'            \item effectLog - The magnitude of the current imbalance
#'            \item vote_dist - A crosstable shopwing the proportion of times that each covariate voted to bias towards control, treatment, or to abstain from voting.
#'     }
#'   Only vote_dist is preserved if \code{simple} is \code{TRUE}
#'
#' @returns an object of class RandomisationScheme
#'
#' @references
#' Johns, Hannah, et al. "Common scale minimal sufficient balance: An improved method for covariate‐adaptive randomization based on the Wilcoxon‐Mann‐Whitney odds ratio statistic." Statistics in Medicine 41.10 (2022): 1846-1861.
#'
#' Italiano, Dominic, et al. "Adaptive Randomization Method to Prevent Extreme Instances of Group Size and Covariate Imbalance in Stroke Trials." Stroke 55.8 (2024): 1962-1972.
csmsb <- function(alpha,
                  biasParameters,
                  voteMethod="pooled",
                  biasMethod="static",
                  weights = NULL,
                  groupImbalanceThreshold=Inf,
                  imbalancePriority = "abstain"
){
  out <- list(
    method = "csmsb",
    settings = list(
      alpha=alpha,
      biasParameters=biasParameters,
      voteMethod = voteMethod,
      biasMethod=biasMethod,
      weights=weights,
      groupImbalanceThreshold=groupImbalanceThreshold,
      imbalancePriority=imbalancePriority
    )
  )

  class(out) <- c("RandomisationScheme",class(out))
  out
}


## CPP Implementations ########################################################
csmsb_allocate_cpp <- function(data,
                               arm,
                               newData,
                               settings,
                               randSequence,
                               offset
                       ){

  # Extract out and process inputs
  rand_number <- nrow(data)+1+offset

  alpha <- settings$alpha
  biasParameters <- settings$biasParameters
  biasMethod <- settings$biasMethod
  voteMethod <- settings$voteMethod
  groupImbalanceThreshold <- settings$groupImbalanceThreshold
  imbalancePriority <- settings$imbalancePriority
  weights <- settings$weights


  if(is.null(alpha)) stop("alpha not provided in settings")

  if(is.null(biasParameters)){
    stop("biasParameters not provided in settings (should be optional if permuted block is used, and warned that it does nothing if supplied in this case)")
  }

  if(is.null(biasMethod)) biasMethod <- "static"
  if(is.null(voteMethod)) voteMethod <- "pooled"
  if(is.null(groupImbalanceThreshold))   groupImbalanceThreshold <- Inf
  if(is.null(imbalancePriority)) imbalancePriority <-  "abstain"

  if(!is.null(weights)){

    if(is.null(names(weights))){
      if(length(weights) == ncol(data)){

        names(weights) <- colnames(data)
        if(length(unique(weights)) != 1) warning("weights were not named and are assumed to correspond to columns in data in order")

      } else {
        stop("length of weights should match columns in data if weights are not named")
      }

    } else {



      if(!is.null(setdiff(colnames(data),names(weights)))){
        stop(
          sprintf("No weight found for data columns %s",
                  paste(setdiff(colnames(data),names(weights)), collapse = ",")
          )
        )
      }

      if(!is.null(setdiff(names(weights),colnames(data)))){
        warning(
          sprintf("Discarding weights %s",
                  paste(setdiff(names(weights),colnames(data)), collapse = ",")
          )
        )
      }

      weights <- weights[colnames(data)]

    }

    if(voteMethod != "pooled") stop("weights are only supported if voteMethod is pooled")
  } else {

    if(is.null(weights)){
      weights <- rep(1,ncol(data))
      names(weights) <- colnames(data)
    }

  }


  if(is.finite(groupImbalanceThreshold)){
    if(biasMethod!="static") stop("groupImbalanceThreshold only supported for CS-MSB (static coin)")
  }




  ## End preamble ##########################

  # Below code is specific to cpp implementation



  # Convert options to magic numbers to pass to C++
  voteMethod <- unname(c(majority=0,weighted=1,pooled=2)[voteMethod])
  if(is.na(voteMethod)) stop("voteMethod not recognised")

  biasMethod <- unname(c(static=0,linear=1)[biasMethod])
  if(is.na(biasMethod)) stop("biasMethod not recognised")

  imbalancePriority <- unname(c(abstain=0,covariate=1,group=2)[imbalancePriority])
  if(is.na(imbalancePriority)) stop("imbalancePriority not recognised")


  # Check that existing arms make sense
  if(length(arm) != nrow(data)) stop("arm must contain randomisations for each row of supplied data")
  if(length(setdiff(arm,0:1))>0) stop("arm must contain values of 0 and 1 only")


  # Check that only one sequence provided
  if("list" %in% class(randSequence$sequence)){
    if("Alloc_seq" %in% class(randSequence$sequence[[1]])){
      stop("Only one allocation sequence should be provided, not a list of them.")
    }
  }



  # C++ randomises the last row of its input
  nObs <- nrow(data)
  data <- rbind(data,newData)

  # Parse data type
  dataType <- nClasses <- rep(NA,ncol(data))
  for(i in 1:ncol(data)){

    if("numeric" %in% class(data[,i])){
      dataType[i] <- 0
      nClasses[i] <- 1
    } else if("factor" %in% class(data[,i]) & "ordered" %in% class(data[,i])){
      dataType[i] <- 1
      nClasses[i] <- 1 # Ordinal data is numeric with extra checks, nClasses should be 1 for this purpose
    } else if("factor" %in% class(data[,i])){
      dataType[i] <- 2
      nClasses[i] <- length(levels(data[,i]))
    } else{
      stop(sprintf("Column %s not numeric, factor or ordered factor",colnames(data)[i]))
    }
  }


  # Generate a C++ friendly representation of the data
  tmpdf <- data
  rownames(tmpdf) <- NULL
  for(i in 1:ncol(data)){
    if("factor" %in% class(data[,i]) & "ordered" %in% class(data[,i])){
      tmpdf[,i] <- as.numeric(tmpdf[,i])
    } else if("factor" %in% class(data[,i])){
      tmpdf[,i] <- as.numeric(tmpdf[,i])-1 # Factors have to be zero indexed
    }
  }
  tmpdf <- as.matrix(tmpdf)
  threshold <- rep(abs(qnorm(p = alpha/2)),ncol(tmpdf))

  # Set arbitrarily high if finite
  if(!is.finite(groupImbalanceThreshold)) groupImbalanceThreshold <- 10

  x <- allocatR:::csmsb_get_decision(nObs = nObs,
                            data = tmpdf,
                            dataType_in = dataType,
                            nClasses_in = nClasses,
                            arm_in = arm,
                            threshold = threshold,
                            parms = biasParameters,
                            voteMethod = voteMethod,
                            biasMethod = biasMethod,
                            weights = weights,
                            groupImbalanceThreshold=groupImbalanceThreshold,
                            imbalancePriority=imbalancePriority
  )

  x$testStatistic <- ifelse(x$testStatistic<0,NA,x$testStatistic)

  if(voteMethod==2){
    # Pooled voting only uses a single test statistic. Trim the excess NAs from output
    x$testStatistic <- x$testStatistic[1]
  }

  x$vote <- ifelse(x$vote<0,NA,x$vote)

  do.call("rbind",
          lapply(1:length(x$imbalance_prospective),function(i){
            out <- c(
              prospective_0=x$imbalance_prospective[[i]][[1]][1],
              prospective_1=x$imbalance_prospective[[i]][[2]][1]
            )
            out <- ifelse(out<=0,NA,out)
            out <- abs(log(out))
            return(out)
          })
  )-> prospectiveImbalance

  effect <- do.call("c",lapply(x$imbalance,function(x){x[1:2]}))
  effect[2*(1:(length(effect)/2))] <- effect[2*(1:(length(effect)/2))]/effect[2*(1:(length(effect)/2))-1]
  effect[2*(1:(length(effect)/2))-1] <- abs(log(effect[2*(1:(length(effect)/2))-1]))

  BD <- x$BD
  seBD <- x$seBD

  details <- list(
    testStatistic = x$testStatistic,
    vote=x$vote,
    imbalance = effect[2*(1:(length(effect)/2))-1],
    seImbalance = effect[2*(1:(length(effect)/2))],
    prospectiveImbalance0 = prospectiveImbalance[,1],
    prospectiveImbalance1 = prospectiveImbalance[,2],
    BD = x$BD,
    seBD = x$seBD,
    probability_covariate=x$probability_covariate,
    probability_group=x$probability_group
  )

  ## Random Assignment ##################################################

  # Get the random number for this allocation

  if("Runif" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence
    roll <- randSequence[rand_number]

  } else if ("Permutedblock" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence
    roll <- randSequence[rand_number,]

    if(roll[1] == 0){

      if(roll[2] == 0 ){

        roll <- (0+(1-biasParameters)) / 2

      } else if (roll[2]==1){

        roll <- ((1-biasParameters) + 0.5) / 2

      } else stop("Something is wrong with randSequence")

    } else if (roll[1] == 1) {

      if(roll[2] == 0 ){

        roll <- ((biasParameters) + 1) / 2

      } else if (roll[2]==1){

        roll <- (0.5 + (biasParameters)) / 2

      } else stop("Something is wrong with randSequence")

    } else {
      stop("Something is wrong with randSequence")
    }
  }

  decision <- ifelse(roll<x$probability,0,1)

  return(list(
    decision = decision,
    prob = x$probability,
    roll = roll,
    biasDirection = x$biasDirection,
    details = details
  ))

}

csmsb_trialSim_cpp <- function(data,
                               settings,
                               randSequence,
                               replay){

  # Extract out and process inputs
  alpha <- settings$alpha
  biasParameters <- settings$biasParameters
  biasMethod <- settings$biasMethod
  voteMethod <- settings$voteMethod
  groupImbalanceThreshold <- settings$groupImbalanceThreshold
  imbalancePriority <- settings$imbalancePriority
  weights <- settings$weights


  if(is.null(alpha)) stop("alpha not provided in settings")

  if(is.null(biasParameters)){
    stop("biasParameters not provided in settings (should be optional if permuted block is used, and warned that it does nothing if supplied in this case)")
  }

  if(is.null(biasMethod)) biasMethod <- "static"
  if(is.null(voteMethod)) voteMethod <- "pooled"
  if(is.null(groupImbalanceThreshold))   groupImbalanceThreshold <- Inf
  if(is.null(imbalancePriority)) imbalancePriority <-  "abstain"


  if(!is.null(weights)){

    if(is.null(names(weights))){
      if(length(weights) == ncol(data)){

        names(weights) <- colnames(data)
        if(length(unique(weights)) != 1) warning("weights were not named and are assumed to correspond to columns in data in order")

      } else {
        stop("length of weights should match columns in data if weights are not named")
      }

    } else {



      if(!is.null(setdiff(colnames(data),names(weights)))){
        stop(
          sprintf("No weight found for data columns %s",
                  paste(setdiff(colnames(data),names(weights)), collapse = ",")
          )
        )
      }

      if(!is.null(setdiff(names(weights),colnames(data)))){
        warning(
          sprintf("Discarding weights %s",
                  paste(setdiff(names(weights),colnames(data)), collapse = ",")
          )
        )
      }

      weights <- weights[colnames(data)]

    }

    if(voteMethod != "pooled") stop("weights are only supported if voteMethod is pooled")
  } else {

    if(is.null(weights)){
      weights <- rep(1,ncol(data))
      names(weights) <- colnames(data)
    }

  }

  if(is.finite(groupImbalanceThreshold)){
    if(biasMethod!="static") stop("groupImbalanceThreshold only supported for CS-MSB (static coin)")
  }

  # Convert options to magic numbers to pass to C++
  voteMethod <- unname(c(majority=0,weighted=1,pooled=2)[voteMethod])
  if(is.na(voteMethod)) stop("voteMethod not recognised")

  biasMethod <- unname(c(static=0,linear=1)[biasMethod])
  if(is.na(biasMethod)) stop("biasMethod not recognised")

  imbalancePriority <- unname(c(abstain=0,covariate=1,group=2)[imbalancePriority])
  if(is.na(imbalancePriority)) stop("imbalancePriority not recognised")

  # Parse data type
  dataType <- nClasses <- rep(NA,ncol(data))
  for(i in 1:ncol(data)){

    if("numeric" %in% class(data[,i])){
      dataType[i] <- 0
      nClasses[i] <- 1
    } else if("factor" %in% class(data[,i]) & "ordered" %in% class(data[,i])){
      dataType[i] <- 1
      nClasses[i] <- 1 # Ordinal data is numeric with extra checks, nClasses should be 1 for this purpose
    } else if("factor" %in% class(data[,i])){
      dataType[i] <- 2
      nClasses[i] <- length(levels(data[,i]))
    } else{
      stop(sprintf("Column %s not numeric, factor or ordered factor",colnames(data)[i]))
    }
  }

  # Generate a C++ friendly representation of the data
  tmpdf <- data
  rownames(tmpdf) <- NULL
  for(i in 1:ncol(data)){
    if("factor" %in% class(data[,i]) & "ordered" %in% class(data[,i])){
      tmpdf[,i] <- as.numeric(tmpdf[,i])
    } else if("factor" %in% class(data[,i])){
      tmpdf[,i] <- as.numeric(tmpdf[,i])-1 # Factors have to be zero indexed
    }
  }
  tmpdf <- as.matrix(tmpdf)

  # Set arbitrarily high if finite
  if(!is.finite(groupImbalanceThreshold)) groupImbalanceThreshold <- 10


  # Generate random numbers for randomisation and feed this list to C++
  if("Runif" %in% class(randSequence$sequence)){

    roll <- randSequence$sequence

  } else if ("Permutedblock" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence

    roll <- rep(NA,nrow(randSequence))
    roll[randSequence[,1] == 0 & randSequence[,2]==0] <- (0+(1-biasParameters)) / 2
    roll[randSequence[,1] == 0 & randSequence[,2]==1] <- ((1-biasParameters) + 0.5) / 2
    roll[randSequence[,1] == 1 & randSequence[,2]==1] <- (0.5 + (biasParameters)) / 2
    roll[randSequence[,1] == 1 & randSequence[,2]==0] <- ((biasParameters) + 1) / 2

    if(any(is.na(roll)))  stop("Something is wrong with randSequence")
  }


  threshold <- rep(abs(qnorm(p = alpha/2)),ncol(tmpdf))

  allocatR:::csmsb_simulate(dataType_in = dataType,nClasses_in = nClasses,
                         data = tmpdf,roll = roll,
                         threshold_in = threshold,
                         parms_in = biasParameters,
                         voteMethod = voteMethod,
                         biasMethod = biasMethod,
                         weights = weights,
                         groupImbalanceThreshold=groupImbalanceThreshold,
                         imbalancePriority=imbalancePriority
  ) -> trial_result

  trial_result$imbalance[trial_result$imbalance==-1] <- NA
  trial_result$vote[trial_result$vote==-1] <- NA
  trial_result$testStatistic[trial_result$testStatistic==-9999] <- NA
  trial_result$testStatistic[is.nan(trial_result$testStatistic)] <- NA

  # Format results

  arm <- factor(trial_result$arm,levels=0:1,labels=c("control","treatment"))

  probabilityLog <- trial_result$randProb
  voteLog <- trial_result$vote
  voteLog[voteLog == -1] <- NA
  colnames(voteLog) <- colnames(data)

  probabilityLog_covariate <- trial_result$randProb_covariate
  probabilityLog_group <- trial_result$randProb_group

  if(voteMethod==2){
    pValLog <- trial_result$testStatistic[,1]
    pValLog[pValLog==-1] <- NA
    pValLog <- 2*pnorm(pValLog,lower.tail = FALSE)
  }else{
    pValLog <- trial_result$testStatistic
    colnames(pValLog) <- colnames(data)
    pValLog <- trial_result$testStatistic
    pValLog[pValLog==-1] <- NA
    pValLog <- 2*pnorm(pValLog,lower.tail = FALSE)
  }


  effectLog <- trial_result$imbalance
  effectLog[effectLog==-1] <- NA
  colnames(effectLog) <- c(t(outer(colnames(data),c("or","se"),"paste",sep="_")))

  # Convert to log scale
  effectLog[,2*(1:(ncol(effectLog)/2))] <- effectLog[,2*(1:(ncol(effectLog)/2))]/effectLog[,2*(1:(ncol(effectLog)/2))-1]
  effectLog[,2*(1:(ncol(effectLog)/2))-1] <- log(effectLog[,2*(1:(ncol(effectLog)/2))-1])

  out <- list(voteLog = voteLog,
              pValLog = pValLog,
              probabilityLog=probabilityLog,
              probabilityLog_covariate=probabilityLog_covariate,
              probabilityLog_group=probabilityLog_group,
              rollLog=roll,
              effectLog= effectLog,
              arm = arm)

  out$imbalance_std <- getImbalance(data,arm,standard=T)
  out$imbalance_genodds <- getImbalance(data,arm,standard = F)
  out$prop_random <- getProportionRandom(probabilityLog)

  out$prop_random_covariate <- getProportionRandom(probabilityLog_covariate)
  out$prop_random_group <- getProportionRandom(probabilityLog_group)

  out$expected_bias <- getExpectedBias(probabilityLog)
  out$vote_dist <- getVoteDist(voteLog)

  return(out)

}


## R Implementations ########################################################
csmsb_allocate_R <- function(data,
                             arm,
                             newData,
                             settings,
                             randSequence,
                             offset){

  # Extract out and process inputs
  rand_number <- nrow(data)+1+offset

  alpha <- settings$alpha
  biasParameters <- settings$biasParameters
  biasMethod <- settings$biasMethod
  voteMethod <- settings$voteMethod
  groupImbalanceThreshold <- settings$groupImbalanceThreshold
  imbalancePriority <- settings$imbalancePriority
  weights <- settings$weights


  if(is.null(alpha)) stop("alpha not provided in settings")

  if(is.null(biasParameters)){
    stop("biasParameters not provided in settings (should be optional if permuted block is used, and warned that it does nothing if supplied in this case)")
  }

  if(is.null(biasMethod)) biasMethod <- "static"
  if(is.null(voteMethod)) voteMethod <- "pooled"
  if(is.null(groupImbalanceThreshold))   groupImbalanceThreshold <- Inf
  if(is.null(imbalancePriority)) imbalancePriority <-  "abstain"

  if(!is.null(weights)){

    if(is.null(names(weights))){
      if(length(weights) == ncol(data)){

        names(weights) <- colnames(data)
        if(length(unique(weights)) != 1) warning("weights were not named and are assumed to correspond to columns in data in order")

      } else {
        stop("length of weights should match columns in data if weights are not named")
      }

    } else {



      if(!is.null(setdiff(colnames(data),names(weights)))){
        stop(
          sprintf("No weight found for data columns %s",
                  paste(setdiff(colnames(data),names(weights)), collapse = ",")
          )
        )
      }

      if(!is.null(setdiff(names(weights),colnames(data)))){
        warning(
          sprintf("Discarding weights %s",
                  paste(setdiff(names(weights),colnames(data)), collapse = ",")
          )
        )
      }

      weights <- weights[colnames(data)]

    }

    if(voteMethod != "pooled") stop("weights are only supported if voteMethod is pooled")
  } else {

    if(is.null(weights)){
      weights <- rep(1,ncol(data))
      names(weights) <- colnames(data)
    }

  }

  if(is.finite(groupImbalanceThreshold)){
    if(biasMethod!="static") stop("groupImbalanceThreshold only supported for CS-MSB (static coin)")
  }

  ## End preamble #######################################################

  # Initialise outputs

  decision <- NA
  biasDirection <- NA

  votes <- rep(NA,ncol(newData))
  pVals <- rep(NA,ncol(newData))
  effect <- rep(NA,ncol(newData))
  effectSe <- rep(NA,ncol(newData))
  effectLower <- rep(NA,ncol(newData))
  effectUpper <- rep(NA,ncol(newData))
  prospective_effect0 <- rep(NA,ncol(newData))
  prospective_effect1 <- rep(NA,ncol(newData))

  pooled_L <- NA
  se_L <- NA
  pool_p <- NA

  for( j in 1:length(votes)){

    if(("factor" %in% class(newData[,j])) & !("ordered" %in% class(newData[,j])))
    {
      # If numeric, run genodds on binary "does covariate match new variable"
      if(alpha < 0)
      {
        # If we've fed in a nonsensical alpha, suppress any warnings as a result of this
        fit <- tryCatch(suppressWarnings({
          genodds::genodds(data[,j]==newData[,j],arm,alpha = alpha)
        }),
        error=function(e){list(pooled_p=NA)})
      }
      else
      {
        fit <- tryCatch(genodds::genodds(data[,j]==newData[,j],arm,alpha = alpha),error=function(e){list(pooled_p=NA)})
      }

      pVals[j] <- fit$pooled_p

      if(!is.na(fit$pooled_p))
      {

        prospective_fit1 <- genodds::genodds(c(unlist(data[,j]==newData[,j]),unlist(newData[,j]==newData[,j])),c(arm,1))
        prospective_fit0 <- genodds::genodds(c(unlist(data[,j]==newData[,j]),unlist(newData[,j]==newData[,j])),c(arm,0))

        prospective_effect1[j] <- prospective_fit1$pooled_lnodds
        prospective_effect0[j] <- prospective_fit0$pooled_lnodds

        if(abs(prospective_fit1$pooled_lnodds) < abs(prospective_fit0$pooled_lnodds))
        {
          # Assigning to 1 pushes us closer to null
          votes[j] <- 1
        }
        else if (abs(prospective_fit1$pooled_lnodds) > abs(prospective_fit0$pooled_lnodds))
        {
          # Assigning to 1 pushes us away from null
          votes[j] <- 0
        }
        else
        {
          votes[j] <- NA
        }

        if(!is.na(votes[j]))
        {
          effect[j] <- (-1)^(votes[j]+1) * abs(fit$pooled_lnodds)
          effectSe[j] <- fit$pooled_SElnodds
          effectLower[j] <- fit$pooled_lnconf.int[1]
          effectUpper[j] <- fit$pooled_lnconf.int[2]
        }

      }
    }
    else
    {
      if(alpha < 0)
      {
        # If we've fed in a nonsensical alpha, suppress any warnings as a result of this
        fit <- tryCatch(suppressWarnings({
          genodds::genodds(unlist(data[,j]),arm,alpha = alpha)
        }),
        error=function(e){list(pooled_p=NA)})
      }
      else
      {
        fit <- tryCatch(genodds::genodds(unlist(data[,j]),arm,alpha = alpha),error=function(e){list(pooled_p=NA)})
      }

      pVals[j] <- fit$pooled_p

      if(!is.na(fit$pooled_p))
      {

        #Prospectively look at the imbalance based on this decision
        # If we assign the new record to 1 does it make the odds ratio larger or smaller?

        prospective_fit1 <- genodds::genodds(c(unlist(data[,j]),unlist(newData[,j])),c(arm,1))
        prospective_fit0 <- genodds::genodds(c(unlist(data[,j]),unlist(newData[,j])),c(arm,0))

        prospective_effect1[j] <- prospective_fit1$pooled_lnodds
        prospective_effect0[j] <- prospective_fit0$pooled_lnodds

        if(abs(prospective_fit1$pooled_lnodds) < abs(prospective_fit0$pooled_lnodds))
        {
          # Assigning to 1 pushes us closer to null
          votes[j] <- 1
        }
        else if (abs(prospective_fit1$pooled_lnodds) > abs(prospective_fit0$pooled_lnodds))
        {
          # Assigning to 1 pushes us away from null
          votes[j] <- 0
        }
        else
        {
          votes[j] <- NA
        }

        if(!is.na(votes[j]))
        {
          effect[j] <- ifelse(is.na(votes[j]),0,(-1)^(votes[j]+1) * abs(fit$pooled_lnodds))
          effectSe[j] <- fit$pooled_SElnodds
          effectLower[j] <- fit$pooled_lnconf.int[1]
          effectUpper[j] <- fit$pooled_lnconf.int[2]
        }
        else
        {
          effect[j] <- 0
          effectSe[j] <- fit$pooled_SElnodds
          effectLower[j] <- fit$pooled_lnconf.int[1]
          effectUpper[j] <- fit$pooled_lnconf.int[2]
        }
      }
    }
  }

  # Step 2: combine covariates and determine coin bias
  if(voteMethod %in% c("majority","weighted")){
    # Count significant votes and get probability
    if(voteMethod=="majority")
    {
      # Simple voting
      sum0 <- sum(votes==0 & pVals<=alpha ,na.rm=T)
      sum1 <- sum(votes==1 & pVals<=alpha ,na.rm=T)
    }
    else if(voteMethod=="weighted")
    {
      #Weighted voting
      sum0 <- sum(abs(effect)[which(votes==0 & pVals<=alpha )],na.rm=T)
      sum1 <- sum(abs(effect)[which(votes==1 & pVals<=alpha )],na.rm=T)
    }

    if(sum0>sum1){
      prob <- biasParameters[1]
      biasDirection <- 0
    } else if(sum1 > sum0){
      prob <- 1-biasParameters[1]
      biasDirection <- 1
    } else {
      prob <- 0.5
      biasDirection <- NA
    }

  } else if (voteMethod=="pooled"){

    # Pooled voting using inverse variance weighting
    # ln(odds) < 0 indicates a vote for 0, ln(odds)>0 indicates a vote for 1

    if(is.null(weights))
    {
      weights <- rep(1,length(effect))
    }

    pooled_L <- sum(weights * effect / effectSe^2 ,na.rm = FALSE)/sum(weights/effectSe^2,na.rm = FALSE)
    se_L <- sum(weights^2/effectSe^2)/(sum(weights/effectSe^2)^2)

    pool_p <- pnorm(abs(pooled_L),mean = 0,sd=sqrt(se_L),lower.tail = FALSE)*2

    if(is.na(pool_p)) pool_p <- Inf

    if(pool_p <= alpha)
    {
      # Get probability of voting 0.
      # scaleFactor controls how extreme the biasing is

      if(biasMethod=="static")
      {
        if(pooled_L<0)
        {
          prob <- biasParameters[1]
        }
        else if (pooled_L>0)
        {
          prob <- 1-biasParameters[1]
        }
        else
        {
          prob <- 0.5
        }

      }
      else if(biasMethod=="linear")
      {
        prob <- 0.5+ biasParameters[1]*abs(pooled_L)
      }
      else
      {
        stop("biasMethod not recognised")
      }

      biasDirection <- ifelse(prob>0.5,0,1)


    }
    else
    {
      prob <- 0.5
      biasDirection <- NA
    }
  }


  # Step 3: Estimate group size control and get coin strength for this outcome too
  # Group imbalance is built with a fixed coin in mind

  if(length(arm) > 0){
    p0 <- sum(arm==0)/length(arm)
    p1 <- sum(arm==1)/length(arm)
  } else {
    p0 <- 0
    p1 <- 0
  }

  groupImbalance <- abs(p0-p1)
  groupImbalanceProb <- 0.5

  if(groupImbalance > groupImbalanceThreshold){
    if(p0>p1) groupImbalanceProb <- 1-biasParameters[1]
    if(p0<p1) groupImbalanceProb <- biasParameters[1]
  }

  covBalanceProb <- prob

  if(covBalanceProb!=groupImbalanceProb){

    if(covBalanceProb==0.5){ # covImbalance doesn't care

      prob <- groupImbalanceProb

    } else {
      if(groupImbalanceProb==0.5){ # groupImbalance doesn't care
        prob <- covBalanceProb
      } else {
        # Both covariate balance and group imbalance care
        if(imbalancePriority == "abstain"){
          prob <- 0.5
        } else if(imbalancePriority == "covariate"){
          prob <- covBalanceProb
        } else if(imbalancePriority == "group"){
          prob <- groupImbalanceProb
        }

      }
    }
  }

  details <- list(
    pVals = pVals,
    vote = votes,
    imbalance = effect,
    seImbalance = effectSe,
    prospectiveImbalance0=prospective_effect0,
    prospectiveImbalance1=prospective_effect1,
    BD = pooled_L,
    seBD = se_L,
    BD_pVal = pool_p,
    covBalanceProb=covBalanceProb,
    groupImbalanceProb=groupImbalanceProb
  )


  ## Random Assignment ##################################################

  # Get the random number for this allocation

  if("Runif" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence
    roll <- randSequence[rand_number]

  } else if ("Permutedblock" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence
    roll <- randSequence[rand_number,]

    if(roll[1] == 0){

      if(roll[2] == 0 ){

        roll <- (0+(1-biasParameters)) / 2

      } else if (roll[2]==1){

        roll <- ((1-biasParameters) + 0.5) / 2

      } else stop("Something is wrong with randSequence")

    } else if (roll[1] == 1) {

      if(roll[2] == 0 ){

        roll <- ((biasParameters) + 1) / 2

      } else if (roll[2]==1){

        roll <- (0.5 + (biasParameters)) / 2

      } else stop("Something is wrong with randSequence")

    } else {
      stop("Something is wrong with randSequence")
    }
  }

  decision <- ifelse(roll<prob,0,1)

  return(
    list(
    decision = decision,
    prob = prob,
    roll = roll,
    biasDirection = biasDirection,
    details = details
  ))

}
