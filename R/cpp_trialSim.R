#' Trial Simulation
#'
#' @description Simulate a trial using the MSB, CS-MSB or CSSize-MSB algorithms
#'
#' @usage
#' trialSim(data, outcomes=NULL, primaryAnalysis = NULL,
#'          n=NULL, replay=FALSE,
#'          alpha=0.05,
#'          method='CSMSB',
#'          voteMethod = NULL,
#'          biasMethod = "static",
#'          weights=NULL,
#'          biasParameters,
#'          groupImbalanceThreshold=Inf,
#'          imbalancePriority="abstain",
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
#' @param trng an optional parallel random number generator from the rTRNG package. If an integer is supplied instead, a yarn3 engine will be generated using this integer as the seed.
#' @param jump_trng a nonnegative integer specifying how many steps the trng engine should be jumped forward before it is used.
#'
#' @details
#'
#' See \code{?decideBias} for algorithm details.
#'
#' User generated function for analysis should take the form \code{primaryAnalysis(cov_df,arm,outcome)}, where \code{cov_df} contains a data frame of covariates,\
#' \code{arm} is a vector containing treatment assignments, and \code{outcome} is a vector containing outcome data.
#' While this function can return any value when run alone, if it is to be run through \code{trialOC()} it should return a named vector of values.
#'
#' If \code{replay=TRUE}, then \code{n} is replaced with \code{nrow(data)}.
#'
#' \code{outcomes} must contain a named list which contains two \code{n}\times\code{ncol(data)} matrices.
#' \code{outcomes[[i]][j][k]} contains the probability of patient \code{i} having outcome \code{k} under treatment \code{i}.
#' The names of the list must be \code{"control"} and \code{"treat"}. This means that at present, only binary, ordinal and categorical outcomes are supported.
#' @return  A list containing the following:
#' \itemize{
#'     \item voteLog - A \code{x}\times\code{ncol(data)} matrix containing the votes (0,1,NA) for each covariate at each randomisation. See \code{?decideBias} for details.
#'     \item pValLog - A \code{x}\times\code{ncol(data)} matrix containing the p-values for each covariate vote at each randomisation. See \code{?decideBias} for details.
#'     \item probabilityLog - A vector containing the probability of assignment to control arm at each randomisation
#'     \item rollLog - A vector containing the random uniform number used to determine treatment assignment at each randomisation
#'     \item effectLog -  A \code{x}\times\code{ncol(data)} matrix containing the imbalance for each covariate at each randomisation. See \code{?decideBias} for details.
#'     \item sim_df - A \code{x}\times\code{ncol(data)} data frame containing the covariate data from the simulation
#'     \item arm - A numeric vector containing arm assignments for the simulation (0=control, 1=treatment)
#'     \item outcome  - A vector containing the outcomes for the simulation
#'     \item analysisResults  - An object containing the results of running \code{primaryAnalysis(sim_df,arm,outcome)}
#'     \item alpha  - The supplied \code{alpha}. Retained for other functions (likely depreciated)
#'     \item levels - The levels for each covariate (NA if not factor). Retained for other functions (likely depreciated)
#'     \item imbalance_std - A named vector giving the imbalance in each covariate at the end of the simulation
#'     \item imbalance_genodds - A named vector giving the imbalance in each covariate at the end of the simulation. Nominal covariates are binarised for each value.
#'     \item prop_random - The proportion of purely random assignments
#'     \item expected_bias - The expected deviation from purely random assignment
#'     \item vote_dist - A matrix giving the proportion of votes for each covariate
#' }

trialSim <- function(data,
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

  # Internal Functions for summarising results ################

  getImbalance <- function(covariates,arm,standard=TRUE){
    if(standard)
    {
      statistic <- rep(NA,ncol(covariates))
      pvals <- rep(NA,ncol(covariates))
      for(j in 1:ncol(covariates))
      {

        if("numeric" %in% class(covariates[,j]))
        {
          # Use cohen's d to measure imbalance (pool SD)
          statistic[j] <- (mean(covariates[which(arm=="treatment"),j]) -  mean(covariates[which(arm=="control"),j]))/
            sqrt( (
              sd(covariates[which(arm=="treatment"),j])*(sum(arm=="treatment")-1) +
                sd(covariates[which(arm=="control"),j])*(sum(arm=="control")-1)
            )/(
              nrow(covariates)-1
            )
            )

          pvals[j] <- t.test(covariates[which(arm=="treatment"),j], covariates[which(arm=="control"),j])$p.value

        }
        else if(("factor" %in% class(covariates[,j]))) # We run chisq for nominal and ordinal, because that's what we can do with base MSB
        {
          statistic[j] <- suppressWarnings({chisq.test(table(arm,covariates[,j]))$statistic})
          pvals[j] <- suppressWarnings({chisq.test(table(arm,covariates[,j]))$p.value})
        }
        # else if(("ordered" %in% class(covariates[,j])))
        # {
        #
        #   #wilcox test wants numeric
        #   vals <- as.numeric(covariates[,j])
        #
        #   statistic[j] <- unname(wilcox.test(vals[arm=="treatment"], vals[arm=="control"])$statistic)/(sum(arm=="treatment") * sum(arm=="control"))
        #   pvals[j] <- wilcox.test(vals[arm=="treatment"],vals[arm=="control"])$p.value
        # }
      }
      names(statistic) <- paste("imbalance_",colnames(covariates), sep="")
      names(pvals) <- paste("p_",colnames(covariates), sep="")
    }
    else
    {
      # The size of output is not fixed for genodds so build it dynamically
      statistic <- {}
      pvals <- {}

      for(j in 1:ncol(covariates))
      {
        if(("numeric" %in% class(covariates[,j])) | ("ordered" %in% class(covariates[,j]) ))
        {
          statistic <- c(statistic,genodds::genodds(covariates[,j],arm)$pooled_lnodds)
          pvals <- c(pvals,genodds::genodds(covariates[,j],arm)$pooled_p)
          names(statistic)[length(statistic)] <- colnames(covariates)[j]
          names(pvals)[length(statistic)] <- colnames(covariates)[j]
        }
        else if(("factor" %in% class(covariates[,j])) & !("ordered" %in% class(covariates[,j])))
        {
          for(k in levels(covariates[,j]))
          {
            statistic <- suppressWarnings(c(statistic,genodds::genodds(covariates[,j]==k,arm)$pooled_lnodds))
            pvals <- suppressWarnings(c(pvals,genodds::genodds(covariates[,j]==k,arm)$pooled_p))

            names(statistic)[length(statistic)] <- paste(colnames(covariates)[j],k,sep=":")
            names(pvals)[length(pvals)] <- paste(colnames(covariates)[j],k,sep=":")
          }
        }
      }

      names(statistic) <- paste("imbalance_",names(statistic),sep="")
      names(pvals) <- paste("p_",names(pvals),sep="")
    }
    return(c(statistic,pvals))
  }


  getProportionRandom <- function(x) {
    mean(abs(x==0.5))
  }

  getExpectedBias <- function(x) {
    mean(abs(x-0.5))
  }

  getVoteDist <- function(voteLog) {

    sapply(colnames(voteLog),function(i){
      c(vote0 = sum(voteLog[,i]==0,na.rm = T),
        vote1 = sum(voteLog[,i]==1,na.rm = T),
        abstain = sum(is.na(voteLog[,i]))
      )/nrow(voteLog)
    })
  }


  # Set up random number generator if one was fed through

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

  # Input processing ###########################################################

  if(is.finite(groupImbalanceThreshold))
  {
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


  # Perform resampling if needed

  rownames(data) <- NULL
  if(replay){
    if(!is.null(n)) warning("Input n is ignored when replay is set to TRUE")
    n <- nrow(data)
  } else {
    if(is.null(n)) stop("Must specify sample size n")

    # If we're not replaying the dataset, we need to sample rows to feed through
    if(is.null(trng)){
      rows <- ceiling(runif(n,min=0,max=nrow(data)))
    } else {
      rows <- ceiling(rTRNG::runif_trng(n,min=0,max=nrow(data),engine = trng))
    }

    data <- data[rows,]

    # If we supplied outcomes, we need to do the same for outcome generation
    if(!is.null(outcomes)){
      outcomes <- lapply(outcomes,function(x,rows){
        x[rows,]
      },rows=rows)
    }
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

  if(!is.null(outcomes)){
    # if(is.null(outcomeType)) stop("Must specify type of outcome (numeric, factor, logical, or resample")
    if(!("control" %in% names(outcomes) & "treatment" %in% names(outcomes))) stop("outcomes list must contain named elements 'control' and 'treatment'")
    if(length(setdiff(names(outcomes),c("control","treatment")))) warning("outcomes other than 'control' and 'treatment' ignored")
  }


  #############################################################################

  # At present, the only valid approach is CSMSB
  # Simple randomisation and MSB will be added later

  # Convert options to magic numbers to pass to C++

  # Process inputs and prepare them for dropping to C++

  # data,
  # outcomes=NULL,
  # outcomeType = NULL,
  # primaryAnalysis = NULL,
  # n=100,alpha=0.05,
  # method='CSMSB',
  # voteMethod = 0,
  # fixedBias = TRUE,
  # xi=0.6,
  # scaleFactor=1,
  # replay=FALSE

  # data <- mtcars[,c("mpg","am","gear")]
  # data$gear <- factor(df$gear)


  # Generate random numbers for randomisation and feed this list to C++
  if(is.null(trng)){
    roll <- runif(n)
  } else {
    roll <- rTRNG::runif_trng(n,engine = trng)
  }




  if(method =="MSB"){

    # threshold depends on the method
    threshold <- do.call("rbind",lapply(1:nrow(tmpdf),function(rand_i){
      sapply(1:ncol(tmpdf), function(i, alpha, rand_i){

        if(dataType[i] == 0){
          if(rand_i > 3){
            return(abs(qt(1-alpha/2,df = rand_i-1-2)))
          } else {
            return(Inf)
          }
        }

        if(dataType[i] == 2) return(qchisq(1-alpha,df = (length(unique(tmpdf[,i]))-1)))

      },
      alpha=alpha,rand_i=rand_i)
    }))


    csmsb:::simulate_msb(dataType_in = dataType,nClasses_in = nClasses,
                         data = tmpdf,roll = roll,
                         threshold_in = threshold,
                         parms = biasParameters) -> trial_result

    trial_result$vote[trial_result$vote==-1] <- NA
    trial_result$testStatistic[trial_result$testStatistic==-9999] <- NA
    trial_result$testStatistic[is.nan(trial_result$testStatistic)] <- NA

  } else if(method == "CSMSB"){

    threshold <- rep(abs(qnorm(p = alpha/2)),ncol(tmpdf))

    csmsb:::simulate_csmsb(dataType_in = dataType,nClasses_in = nClasses,
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

  } else {
    stop("Method not supported! Must be either MSB or CSMSB")
  }

  # browser()

  # Format results

  arm <- factor(trial_result$arm,levels=0:1,labels=c("control","treatment"))

  probabilityLog <- trial_result$randProb
  voteLog <- trial_result$vote
  voteLog[voteLog == -1] <- NA
  colnames(voteLog) <- colnames(data)

  if(method == "CSMSB"){
    probabilityLog_covariate <- trial_result$randProb_covariate
    probabilityLog_group <- trial_result$randProb_group
  } else {
    probabilityLog_covariate <- NULL
    probabilityLog_group <- NULL
  }

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

  # effectLog only exists for CSMSB

  if(method == "CSMSB"){
    effectLog <- trial_result$imbalance
    effectLog[effectLog==-1] <- NA
    colnames(effectLog) <- c(t(outer(colnames(data),c("or","se"),"paste",sep="_")))

    # Convert to log scale
    effectLog[,2*(1:(ncol(effectLog)/2))] <- effectLog[,2*(1:(ncol(effectLog)/2))]/effectLog[,2*(1:(ncol(effectLog)/2))-1]
    effectLog[,2*(1:(ncol(effectLog)/2))-1] <- log(effectLog[,2*(1:(ncol(effectLog)/2))-1])

  } else {
    effectLog <- NULL
  }



  # We need metadata if we batch run this so grab it now before we return


  # I literally can't remember what this does but it looks important so leave it for now
  allLevels <- lapply(1:ncol(data),function(x){if("factor" %in% class(data[,x])) return(levels(data[,x])) else return(NA)})

  if(!is.null(outcomes)){
    if(!is.null(trng)) stop("Simulation of outcomes not currently supported by trng. Use R's default rng engine")
    outcome <- sapply(1:n,function(i,outcomes,arm){
      out <- sample(colnames(outcomes[[arm[i]]]),1,prob = outcomes[[arm[i]]][i,])
      return(out)
    },outcomes=outcomes,arm=arm)
  } else {
    outcome <- NULL
  }


  out <- list(voteLog = voteLog,
              pValLog = pValLog,
              probabilityLog=probabilityLog,
              probabilityLog_covariate=probabilityLog_covariate,
              probabilityLog_group=probabilityLog_group,
              rollLog=roll,
              effectLog= effectLog,
              sim_df = data,
              arm = arm,
              outcome=outcome,
              alpha=alpha,
              levels=allLevels)

  out$imbalance_std <- getImbalance(data,arm,standard=T)
  out$imbalance_genodds <- getImbalance(data,arm,standard = F)
  out$prop_random <- getProportionRandom(probabilityLog)

  if(method=="CSMSB"){
    out$prop_random_covariate <- getProportionRandom(probabilityLog_covariate)
    out$prop_random_group <- getProportionRandom(probabilityLog_group)
  }

  out$expected_bias <- getExpectedBias(probabilityLog)
  out$vote_dist <- getVoteDist(voteLog)

  # This is bad and dumb but treating this as a factor is a bad idea for
  # the purposes of the primaryAnalysis function. Let the user handle it.
  arm <- as.numeric(arm)-1

  if(!is.null(primaryAnalysis)){
    analysisResults <- primaryAnalysis(data,arm,outcome)
  } else {
    analysisResults <- NULL
  }

  out$analysisResults <- analysisResults

  return(out)
}
