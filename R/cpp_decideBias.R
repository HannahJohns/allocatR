#' Randomise a participant
#'
#' @description Randomise a participant according to the MSB, CS-MSB or CSSize-MSB algorithms.
#'
#' @usage decideBias(data, arm, newRecord,
#'                   method=c('CSMSB', 'MSB'),
#'                   voteMethod=NULL,
#'                   biasMethod=c("static", "linear"),
#'                   weights = NULL,
#'                   alpha,
#'                   biasParameters,
#'                   groupImbalanceThreshold=Inf,
#'                   imbalancePriority=c("abstain", "covariate", "group")
#'                   trng = NULL,
#'                   jump_trng = 0
#'                   )
#'
#' @param data A data frame of covariates for current records. 0 rows if first randomisation.
#' @param arm  A vector of arm assignments for current records. Length 0 if first randomisation.
#' @param newRecord A 1 by \code{ncol(data)} data frame of covariates for the new record.
#' @param method A character string indicating the type of method to use.  Defaults to CSMSB.
#' @param voteMethod A character string indicating how votes should be combined. Defaults to "majority" for MSB and "pooled" for CSMSB.
#' @param biasMethod A character string indicating how a bias coin probability should be generated. Defaults to "static"
#' @param weights A vector of length \code{ncol(data)} giving the weights for each column in order. Defaults to even weights if not supplied.
#' @param alpha A number indicating the level of significance for when MSB/CSMSB intervenes.
#' @param biasParameters A numeric vector providing parameters to the indicated bias method
#' @param groupImbalanceThreshold Group allocation imbalance threshold required to trigger CSSize-MSB allocation rules. By default will never trigger.
#' @param imbalancePriority Priority rule for CSSize-MSB. Defaults to "abstain"
#' @param trng an optional parallel random number generator from the rTRNG package. If an integer is supplied instead, a yarn3 engine will be generated using this integer as the seed.
#' @param jump_trng a nonnegative integer specifying how many steps the trng engine should be jumped forward before it is used.
#'
#' @details
#'
#' Variable type is read from the class of each column in \code{data}.
#' Factors are assumed categorical, ordered factors are assumed ordinal
#' and numeric are assumed continuous.
#'
#' When \code{method=="MSB"}, the Minimal Sufficient Balance algorithm is used.
#' As implemented here, it uses the Chi-square test to evaluate imbalance on
#' categorical variables, and the two group z-test (assuming equal variance) to
#' evaluate imbalance on continuous variables. This is done even when the
#' chi-square test may be inaccurate due to low cell counts or where the sample variance is large.
#' Ordinal variables are not supported by this algorithm.
#'
#' If the MSB algorithm is used, \code{voteMethod} must be \code{"majority"} and
#' \code{biasMethod} must be \code{"static"}
#'
#' When \code{method=="CSMSB"}, the Common Scale Minimal Sufficient Balance
#' algorithm is used. It always uses the Wilcoxon-Mann-Whitney odds ratio to
#' evaluate imbalance. For categorical variables, it uses a binary representation
#' of the variable, which is 1 if it matches the value in \code{newRecord} or 0 otherwise.
#'
#' \code{biasParameters} should be either a single number or a vector of numbers
#' depending on \code{biasMethod}
#'
#' \enumerate{
#' \item static - a single number between 0.5 and 1 giving the strength of the biased coin
#' \item linear - a single number between 1 and infinity that scales the imbalance measure to calculate a biased coin strength according to the formula 0.5+imbalance*biasParameters
#' }
#'
#' When using majority or weighted voting, there are a few quirks to this implementation.
#' \enumerate{
#'     \item When checking each variable, we assume that a vote is always valid and log the vote as if it is.
#'     \item We remove invalid votes based on \code{alpha}, immediately before we tally votes.
#' }
#' This is mainly done so we can log small imbalances for understanding trial behaviour over time.
#'
#' The \code{details} list returned by this function varies in contents depending on the method used.
#'
#' If \code{method=="MSB"} then it will contain the following:
#'
#' \itemize{
#'  \item testStatistic - the test statistics used for each
#'  \item vote - the votes placed by each variable
#' }
#'
#' If \code{method=="CSMSB"} then it will contain the following:
#'
#' \itemize{
#'  \item testStatistic - The test statistic(s) used by the CSMSB algorithm
#'  \item vote - The votes placed by each variable
#'  \item imbalance - The current imbalance on each variable
#'  \item seImbalance - The standard error of the current imbalance on each variable
#'  \item prospectiveImbalance0 - The imalance on each variable if the new participant is assigned to control
#'  \item prospectiveImbalance1 - The imalance on each variable if the new participant is assigned to intervention
#'  \item BD - The pooled measure of imbalance
#'  \item seBD - The standard error of the pooled measure of imbalance
#'  }
#'
#' @return  A list containing the following:
#' \itemize{
#'     \item decision - A randomisation decision (0 or 1)
#'     \item prob - The probability of assignment to arm 0
#'     \item biasDirection - The direction of the assignment bias (0 or 1)
#'     \item details - A list containing further information about the randomisation decision. The specific contents of this list depend on the method. See Details above.
#'}
decideBias <- function(data,
                       arm,
                       newRecord,
                       method='CSMSB',
                       voteMethod=NULL,
                       biasMethod="static",
                       weights = NULL,
                       alpha,
                       biasParameters,
                       groupImbalanceThreshold=Inf,
                       imbalancePriority="abstain",
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

  if(is.finite(groupImbalanceThreshold)){
    # We can probably do away with the first two of these restrictions but the method was proposed in this context
    if(method != "CSMSB") stop("groupImbalanceThreshold only supported for CS-MSB (pooled imbalance with static coin)")
    if(voteMethod!="pooled") stop("groupImbalanceThreshold only supported for CS-MSB (pooled imbalance with static coin)")
    if(biasMethod!="static") stop("groupImbalanceThreshold only supported for CS-MSB (pooled imbalance with static coin)")
  }

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

  # C++ randomises the last row of its input
  nObs <- nrow(data)
  data <- rbind(data,newRecord)

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



  allowed_engines <- c("Rcpp_yarn2",
                       "Rcpp_yarn3",
                       "Rcpp_yarn4",
                       "Rcpp_yarn5",
                       "Rcpp_yarn3s",
                       "Rcpp_yarn5s",
                       "Rcpp_lcg64",
                       "Rcpp_lcg64_shift",
                       "Rcpp_mrg2",
                       "Rcpp_mrg3",
                       "Rcpp_mrg4",
                       "Rcpp_mrg5",
                       "Rcpp_mrg3s",
                       "Rcpp_mrg5s"
  )

  if(!is.null(trng)){

    if(!("rTRNG" %in% rownames(installed.packages()))){
      stop("rTRNG package must be installed to use trng option")
    }

    # Allow passing of a seed, instead of a fully specified engine
    if(is.numeric(trng)) {

      if(as.integer(trng) != trng){
        stop("trng must be a whole number if a seed is being specified")
      }

      trng <- rTRNG::yarn3$new(trng)

    }

    if(!(class(trng) %in% allowed_engines)){

      stop(sprintf("%s is not an allowed rTRNG parallel engine. trng must be an integer, or any of the following types: %s",
                   class(trng),
                   paste(gsub("Rcpp_","",allowed_engines), collapse=", ")
      ))

    }


    # Jump seed state forward. Use multiple jumps if neccesary.
    while(jump_trng > .Machine$integer.max){
      trng$jump(.Machine$integer.max)
      jump_trng <- jump_trng - .Machine$integer.max
    }

    if(jump_trng > 0){
      trng$jump(jump_trng)
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

  if(method =="MSB"){

    # threshold depends on the method
    threshold <- sapply(1:ncol(tmpdf), function(i, alpha){
      if(dataType[i] == 0) return(abs(qt(1-alpha/2,df = nrow(tmpdf)-1-2)))
      if(dataType[i] == 2) return(qchisq(1-alpha,df = (length(unique(tmpdf[,i]))-1)))
    },
    alpha=alpha)

    x <- csmsb:::msb_get_decision(nObs = nObs,
                              data = tmpdf,
                              dataType_in = dataType,
                              nClasses_in = nClasses,
                              arm_in = arm,
                              threshold = threshold,
                              parms = biasParameters
    )

    x$testStatistic <- ifelse(x$testStatistic<0,NA,x$testStatistic)
    x$vote <- ifelse(x$vote<0,NA,x$vote)

    details <- list(
      testStatistic = x$testStatistic,
      vote=x$vote
    )

  } else if(method == "CSMSB"){

    threshold <- rep(abs(qnorm(p = alpha/2)),ncol(tmpdf))

    # Set arbitrarily high if finite
    if(!is.finite(groupImbalanceThreshold)) groupImbalanceThreshold <- 10

    x <- csmsb:::get_decision(nObs = nObs,
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

  } else {
    stop("Voting method not supported! Must be either MSB or CSMSB")
  }

  if(is.null(trng)){
    roll <- runif(1)
  } else {
    roll <- runif_trng(1,engine=trng)
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
