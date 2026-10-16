#' Minimal Sufficient Balance
#'
#' @description
#' Creates a randomisation scheme for the MSB algorithm
#'
#' @usage
#' msb(alpha,
#'     coinStrength
#' )
#'
#' @param alpha A number indicating the level of significance for when MSB intervenes based on covariate imbalance.
#' @param coinStrength A numeric vector providing the coin bias for MSB. See details.
#'
#' @details
#'
#' When msb() is used with the \link[=allocate]{allocate} function, the details list in its output
#' will contain the following:
#'   \itemize{
#'            \item pVals - p-values indicating imbalance for each covariate
#'            \item vote - The votes placed by each covariate
#'     }
#'
#'
#' When msb() is used with the \link[=trialSim]{trialSim} function, the details list in its output
#' will additionally contain the following:
#'   \itemize{
#'            \item voteLog - the votes placed by each covariate during each allocation
#'            \item pValLog - covariate-level p-value for current imbalance
#'            \item vote_dist - A crosstable shopwing the proportion of times that each covariate voted to bias towards control, treatment, or to abstain from voting.
#'     }
#'   Only vote_dist is preserved if \code{simple} is \code{TRUE}
#'
#' @returns an object of class RandomisationScheme
#'
#' @references
#' Zhao, Wenle, Michael D. Hill, and Yuko Palesch. "Minimal sufficient balance—a new strategy to balance baseline covariates and preserve randomness of treatment allocation." Statistical methods in medical research 24.6 (2015): 989-1002.
#' 
msb <- function(alpha,
                coinStrength
){
  out <- list(
    method = "msb",
    settings = list(
      alpha=alpha,
      coinStrength = coinStrength
    )
  )

  class(out) <- c("RandomisationScheme",class(out))
  out
}

## CPP Implementations ########################################################
msb_allocate_cpp <- function(data,
                             arm,
                             newData,
                             settings,
                             randSequence,
                             offset
                       ){

  # Extract out and process inputs
  rand_number <- nrow(data)+1+offset

  alpha <- settings$alpha
  coinStrength <- settings$coinStrength

  if(is.null(alpha)) stop("alpha not provided in settings")

  if(is.null(coinStrength)){
    stop("coinStrength not provided in settings (should be optional if permuted block is used, and warned that it does nothing if supplied in this case)")
  }

  ## End preamble ##########################

  # Below code is specific to cpp implementation
  
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
      dataType[i] <- 1 # TODO: If we we error later, it's probably because msb doesn't support ordinals and we need to treat it as numeric
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

  x <- allocatR:::CPP_msb_get_decision(nObs = nObs,
                            data = tmpdf,
                            dataType_in = dataType,
                            nClasses_in = nClasses,
                            arm_in = arm,
                            alpha = alpha,
                            parms = coinStrength
  )

  x$pValue <- ifelse(x$pValue>1,NA,x$pValue)
  x$vote <- ifelse(x$vote<0,NA,x$vote)

  details <- list(
    pValue = x$pValue,
    vote=x$vote
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
    details = details
  ))

}

msb_trialSim_cpp <- function(data,
                             settings,
                             randSequence,
                             replay){

  
  # Extract out and process inputs
  alpha <- settings$alpha
  coinStrength <- settings$coinStrength

  if(is.null(alpha)) stop("alpha not provided in settings")

  if(is.null(coinStrength)){
    stop("coinStrength not provided in settings (should be optional if permuted block is used, and warned that it does nothing if supplied in this case)")
  }

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


  allocatR:::CPP_msb_simulate(dataType_in = dataType,nClasses_in = nClasses,
                         data = tmpdf,roll = roll,
                         alpha = alpha,
                         parms = coinStrength
  ) -> trial_result

  trial_result$vote[trial_result$vote==-1] <- NA
  trial_result$pValue[trial_result$pValue>1] <- NA
  trial_result$pValue[is.nan(trial_result$pValue)] <- NA

  # Format results

  arm <- factor(trial_result$arm,levels=0:1,labels=c("control","treatment"))

  probabilityLog <- trial_result$randProb

  voteLog <- trial_result$vote
  voteLog[voteLog == -1] <- NA
  colnames(voteLog) <- colnames(data)

  pValLog <- trial_result$pValue
  colnames(pValLog) <- colnames(data)

  out <- list(voteLog = voteLog,
              pValLog = pValLog,
              probabilityLog=probabilityLog,
              rollLog=roll,
              arm = arm
            )

  out$imbalance_std <- getImbalance(data,arm,standard=T)
  out$imbalance_genodds <- getImbalance(data,arm,standard = F)
  out$prop_random <- getProportionRandom(probabilityLog)

  out$expected_bias <- getExpectedBias(probabilityLog)
  out$vote_dist <- getVoteDist(voteLog)

  return(out)

}


## R Implementations ########################################################
msb_allocate_R <- function(data,
                             arm,
                             newData,
                             settings,
                             randSequence,
                             offset){
            
  # Extract out and process inputs
  rand_number <- nrow(data)+1+offset

  alpha <- settings$alpha
  coinStrength <- settings$coinStrength

  if(is.null(alpha)) stop("alpha not provided in settings")

  if(is.null(coinStrength)){
    stop("coinStrength not provided in settings (should be optional if permuted block is used, and warned that it does nothing if supplied in this case)")
  }

  ## End preamble #######################################################

  # Initialise outputs

  decision <- NA

  votes <- rep(NA,ncol(newData))
  pVals <- rep(NA,ncol(newData))

  n0 <- sum(arm==0,na.rm = T)
  n1 <- sum(arm==1,na.rm = T)

  # Drop empty rows
  if(n0+n1 > 1) {
    data <- data[1:(n0+n1),]
    arm <-  arm[1:(n0+n1)]
  } else {
    data <- data[0,]
    arm <- arm[0]
  }

  for( j in 1:length(votes)){
    if("numeric" %in% class(newData[,j]))
    {

      if(n0 > 1 & n1 > 1){
        m0 <-mean(data[arm==0,j])
        m1 <-mean(data[arm==1,j])
        v0 <- var(data[arm==0,j])
        v1 <- var(data[arm==1,j])
      } else {
        m0 <- m1 <- v0 <- v1 <- NA
      }

      t <- (m0-m1) / sqrt(v1/n1 + v0/n0)

      df <- ( v0/n0 + v1/n1 )^2 /(
        (v0/n0)^2/(n0 -1) + 
        (v1/n1)^2/(n1 -1) 
      )

      pVals[j] <- 2*pt(abs(t),df,lower.tail = FALSE)

      if(!is.na(t)){

        # We're just looking at direction of t here, not significance.
        # All variables get their vote calculated then we just filter for
        # eligible ones at the end.

        if( ( t < 0 &  newData[,j] > m1 ) |  ( t > 0 & newData[,j] < m1 ) ){
          # Vote 0 in either of these situations:

          #         m0----m1----(*)
          #   (*)---m1----m0
          votes[j] <- 0
        } else if( (  t < 0 & newData[,j] < m0 ) |  ( t > 0 & newData[,j] > m0 ) ) {
          # Vote 1 in either of these situations:
          #   (*)----m0----m1
          #          m1----m0----(*)

          votes[j] <- 1
        }

      } # End if is.na(t)
    }# End if class numeric
    else if(("factor" %in% class(newData[,j])))
    {

      # Get crosstab of current counts
      factorLevels <- union(levels(data[,j]),levels(newData[,j]))
      crossTab <- matrix(0,nrow=2,ncol=length(factorLevels))
      for(leveli in 1:length(factorLevels)){
        crossTab[1,leveli] = sum(as.character(data[arm==0,j]) == as.character(factorLevels[leveli]))
        crossTab[2,leveli] = sum(as.character(data[arm==1,j]) == as.character(factorLevels[leveli]))
      }
      colnames(crossTab) <- factorLevels

      expectedTab <- matrix(0,nrow=2,ncol=length(factorLevels))
      for(leveli in 1:length(factorLevels)){
        expectedTab[1,leveli] <- sum(crossTab[,leveli]) * sum(crossTab[1,]) / sum(crossTab)
        expectedTab[2,leveli] <- sum(crossTab[,leveli]) * sum(crossTab[2,]) / sum(crossTab)
      }
      colnames(expectedTab) <- factorLevels

      chisq <- sum((crossTab - expectedTab)^2/expectedTab)
      if(!is.nan(chisq)){

        pVal <- pchisq(q = chisq, df = (length(factorLevels)-1), lower.tail = FALSE)
        pVals[j] <- pVal

        if(crossTab[1,as.character(newData[,j])] < expectedTab[1,as.character(newData[,j])])
        {
          votes[j] <- 0
        } else if(crossTab[2,as.character(newData[,j])] < expectedTab[2,as.character(newData[,j])])
        {
          votes[j] <- 1
        }


      }
      else
      {
        # Set pVal to 2 so that it's never used 
        pVals[j] <- 2
      }

    } # End if categorical

  }# End for column

  
  # Step 2: combine covariates and determine coin bias
  sum0 <- sum(votes==0 & pVals<=alpha ,na.rm=T)
  sum1 <- sum(votes==1 & pVals<=alpha ,na.rm=T)

  if(sum0>sum1){
    prob <- coinStrength
  } else if(sum1 > sum0){
    prob <- 1-coinStrength
  } else {
    prob <- 0.5
  }


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

  # If pVal was >1 then it should be NA
  pVals[pVals>1] <- NA

  details <- list(
    pVals = pVals,
    vote = votes
  )

  return(
    list(
    decision = decision,
    prob = prob,
    roll = roll,
    details = details
  ))

}


