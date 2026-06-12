#' Randomise a patient according to the MSB, CS-MSB or CSSize-MSB algorithms (R implementation)
#'
#' @description Description goes here.
#'
#' @usage R_decideBias(data, arm, newRecord,
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
#' This is a functionally equivalent implementation of \code{decideBias()}, except
#' it is implemented entirely within R. It is therefore much slower in a simulation
#' context, but it does not require the same links to compiled code. It should only be
#' used if decideBias() is returning unresolvable errors (due to e.g. incompatibility with
#' some OS). There is some differnce in output formatting.
#'
#' The \code{details} list returned by this function varies in contents depending on the method used.
#'
#' If \code{method=="MSB"} then it will contain the following:
#'
#' \itemize{
#'  \item pVals - the p-values used for each variable
#'  \item testStatistic - the test statistic used for each variable
#'  \item vote - the votes placed by each variable
#' }
#'
#' If \code{method=="CSMSB"} then it will contain the following:
#' \itemize{
#'  \item pVals - The p values for imbalance on each variable
#'  \item vote - The votes placed by each variable
#'  \item imbalance - The current imbalance on each variable
#'  \item seImbalance - The standard error of the current imbalance on each variable
#'  \item prospectiveImbalance0 - The imalance on each variable if the new participant is assigned to control
#'  \item prospectiveImbalance1 - The imalance on each variable if the new participant is assigned to intervention
#'  \item BD - The pooled measure of imbalance
#'  \item seBD - The standard error of the pooled measure of imbalance
#'  \item BD_pVal = The p-value for the pooled measure of imbalance
#'  }
#'
#' @return  A list containing the following:
#' \itemize{
#'     \item decision - A randomisation decision (0 or 1)
#'     \item prob - The probability of assignment to arm 0
#'     \item biasDirection - The direction of the assignment bias (0 or 1)
#'     \item details - A list containing further information about the randomisation decision. The specific contents of this list depend on the method. See Details above.
#'}
R_decideBias <- function(data,
                         arm,
                         newRecord,
                         method = "CSMSB",
                         voteMethod = NULL,
                         biasMethod = "static",
                         weights=NULL,
                         alpha,
                         biasParameters,
                         groupImbalanceThreshold = Inf,
                         imbalancePriority = "abstain",
                         trng = NULL,
                         jump_trng = 0
                       )
{


  if(method == "CSMSB"){
    if(is.null(voteMethod)) voteMethod <- "pooled"
  } else if(method=="MSB"){
    if(is.null(voteMethod)) voteMethod <- "majority"
    if(voteMethod != "majority") stop("MSB algorithm only supports majority voting")
  }

  # Not all inputs are valid, check these here:


  # method="MSB" should only be used with voteMethod=0 (simple voting)
  if(method == "MSB" & voteMethod !="majority") stop('voteMethod must be "majority" for standard MSB')

  # voteMethod=0 (simple voting) or 1 (weighted voting) should only be used with fixedBias=T
  if(voteMethod %in% c("majority","weighted") & biasMethod != "static") stop('voteMethods of "majority" or "weighted" must have fixedBias=TRUE')


  if(is.finite(groupImbalanceThreshold))
  {
    # We can probably do away with the first two of these restrictions but the method was proposed in this context
    if(method != "CSMSB") stop("groupImbalanceThreshold only supported for CS-MSB (pooled imbalance with static coin)")
    if(voteMethod!="pooled") stop("groupImbalanceThreshold only supported for CS-MSB (pooled imbalance with static coin)")
    if(biasMethod!="static") stop("groupImbalanceThreshold only supported for CS-MSB (pooled imbalance with static coin)")
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


  # Initialise outputs

  decision <- NA
  biasDirection <- NA

  votes <- rep(NA,ncol(newRecord))
  pVals <- rep(NA,ncol(newRecord))
  effect <- rep(NA,ncol(newRecord))
  effectSe <- rep(NA,ncol(newRecord))
  effectLower <- rep(NA,ncol(newRecord))
  effectUpper <- rep(NA,ncol(newRecord))
  prospective_effect0 <- rep(NA,ncol(newRecord))
  prospective_effect1 <- rep(NA,ncol(newRecord))

  pooled_L <- NA
  se_L <- NA
  pool_p <- NA


  # Put checks for input cleanliness here


  # First, get the results
  if(method =="MSB"){

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
      if("numeric" %in% class(newRecord[,j]))
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
        pVals[j] <- 2*pt(abs(t),n0+n1-2,lower.tail = FALSE)

        if(!is.na(t)){

          effect[j] <- t

          # We're just looking at direction of t here, not significance.
          # All variables get their vote calculated then we just filter for
          # eligible ones at the end.

          if( ( t < 0 &  newRecord[,j] > m1 ) |  ( t > 0 & newRecord[,j] < m1 ) ){
            # Vote 0 in either of these situations:

            #         m0----m1----(*)
            #   (*)---m1----m0
            votes[j] <- 0
          } else if( (  t < 0 & newRecord[,j] < m0 ) |  ( t > 0 & newRecord[,j] > m0 ) ) {
            # Vote 1 in either of these situations:
            #   (*)----m0----m1
            #          m1----m0----(*)

            votes[j] <- 1
          }

        } # End if is.na(t)
      }# End if class numeric
      else if(("factor" %in% class(newRecord[,j])))
      {

        # Get crosstab of current counts
        factorLevels <- union(levels(data[,j]),levels(newRecord[,j]))
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
          effect[j] <- chisq

          if(crossTab[1,as.character(newRecord[,j])] < expectedTab[1,as.character(newRecord[,j])])
          {
            votes[j] <- 0
          } else if(crossTab[2,as.character(newRecord[,j])] < expectedTab[2,as.character(newRecord[,j])])
          {
            votes[j] <- 1
          }


        }
        else
        {
          pVals[j] <- 2
        }

      } # End if categorical


    }# End for column

    details <- list(
      vote = votes,
      testStatistic = effect,
      pVals = pVals
    )


  } else if(method=="CSMSB") {

    for( j in 1:length(votes)){

      if(("factor" %in% class(newRecord[,j])) & !("ordered" %in% class(newRecord[,j])))
      {
        # If numeric, run genodds on binary "does covariate match new variable"
        if(alpha < 0)
        {
          # If we've fed in a nonsensical alpha, suppress any warnings as a result of this
          fit <- tryCatch(suppressWarnings({
            genodds::genodds(data[,j]==newRecord[,j],arm,alpha = alpha)
          }),
          error=function(e){list(pooled_p=NA)})
        }
        else
        {
          fit <- tryCatch(genodds::genodds(data[,j]==newRecord[,j],arm,alpha = alpha),error=function(e){list(pooled_p=NA)})
        }

        pVals[j] <- fit$pooled_p

        if(!is.na(fit$pooled_p))
        {

          prospective_fit1 <- genodds::genodds(c(unlist(data[,j]==newRecord[,j]),unlist(newRecord[,j]==newRecord[,j])),c(arm,1))
          prospective_fit0 <- genodds::genodds(c(unlist(data[,j]==newRecord[,j]),unlist(newRecord[,j]==newRecord[,j])),c(arm,0))

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

          prospective_fit1 <- genodds::genodds(c(unlist(data[,j]),unlist(newRecord[,j])),c(arm,1))
          prospective_fit0 <- genodds::genodds(c(unlist(data[,j]),unlist(newRecord[,j])),c(arm,0))

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

  } else {
    stop("Voting method not supported! Must be either MSB or CSMSB")
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
    se_L <- 1/sum(effectSe^-2)

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

  if(is.null(trng)){
    roll <- runif(1)
  } else {
    roll <- runif_trng(1,engine=trng)
  }

  decision <- ifelse(roll<prob,0,1)

  # Return decision and reasoning behind it

  return(list(
    decision = decision,
    prob = prob,
    roll = roll,
    biasDirection = biasDirection,
    details = details
  ))

}




# decideBias(
#   data=mtcars[,c("mpg","cyl","disp")],
#   arm = rep(c(0,1),each=nrow(mtcars)/2),
#   newRecord = mtcars[6,c("mpg","cyl","disp")],
#   method="genodds",
#   voteMethod = 2,
#   alpha=0.34,
#   fixedBias = TRUE,
#   xi=0.6,
#   scaleFactor=1
# )

#
# decideBias(
#   data=mtcars[,c("cyl","mpg","disp")],
#   arm = rep(c(0,1),each=nrow(mtcars)/2),
#   newRecord = mtcars[6,c("cyl","mpg","disp")],
#   method="genodds",
#   voteMethod = 0,
#   alpha=0.34,
#   fixedBias = TRUE,
#   xi=0.6,
#   scaleFactor=1
# )
