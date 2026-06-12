trialSim_minimisation <- function(data, stratum=NULL,
                                  biasStrength,
                                  n=NULL,
                                  replay=FALSE,
                                  balanceGroups = FALSE,
                                  rollBlocks=NULL,
                                  interrims=NULL)
{

  # At present, the only valid approach is CSMSB
  # Simple randomisation and MSB will be added later

  # Convert options to magic numbers to pass to C++
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

          pvals[j] <- tryCatch(t.test(covariates[which(arm=="treatment"),j], covariates[which(arm=="control"),j])$p.value,error=function(e){NA})

        }
        else if(("factor" %in% class(covariates[,j]))) # We run chisq for nominal and ordinal, because that's what we can do with base MSB
        {
          statistic[j] <- tryCatch(chisq.test(table(arm,covariates[,j]))$statistic,error=function(e){NA})
          pvals[j] <- tryCatch(chisq.test(table(arm,covariates[,j]))$p.value,error=function(e){NA})
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
          statistic <- c(statistic,tryCatch(genodds::genodds(covariates[,j],arm)$pooled_lnodds,error=function(e){NA}))
          pvals <- c(pvals,tryCatch(genodds::genodds(covariates[,j],arm)$pooled_p,error=function(e){NA}))
          names(statistic)[length(statistic)] <- colnames(covariates)[j]
          names(pvals)[length(statistic)] <- colnames(covariates)[j]
        }
        else if(("factor" %in% class(covariates[,j])) & !("ordered" %in% class(covariates[,j])))
        {
          for(k in levels(covariates[,j]))
          {
            statistic <- suppressWarnings(c(statistic,tryCatch(genodds::genodds(covariates[,j]==k,arm)$pooled_lnodds,error=function(e){NA})))
            pvals <- suppressWarnings(c(pvals,tryCatch(genodds::genodds(covariates[,j]==k,arm)$pooled_p,error=function(e){NA})))

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

  if(prod(sapply(1:ncol(data),
                 function(i){"factor" %in% class(data[,i])}
                 )
          )==0
  ){
    stop("All covariates must be factors")
  }

  nClasses  <- unname(sapply(colnames(data),function(i){length(unique(data[,i]))}))


  rownames(data) <- NULL
  if(replay){
    if(!is.null(n)) warning("Input n is ignored when replay is set to TRUE")
    n <- nrow(data)

    if(is.null(stratum)){
      stratum <- rep(0,n)
    }
  } else {
    if(is.null(n)) stop("Must specify sample size n")

    # If we're not replaying the dataset, we need to sample rows to feed through
    rows <- sample(1:nrow(data),n,replace = TRUE)
    data <- data[rows,]

    if(is.null(stratum)){
      stratum <- rep(0,n)
    } else {
      stratum <- stratum[rows]
    }
  }

  data_cleaned <- data
  for(i in colnames(data_cleaned)) data_cleaned[,i] <- as.integer(data_cleaned[,i])-1
  data_cleaned <- as.matrix(data_cleaned)

  stratum_cleaned <- as.numeric(factor(stratum))-1


  if(is.null(rollBlocks)){
    roll = runif(nrow(data_cleaned))
  } else {

    if(length(rollBlocks)!= n) stop("rollBlocks must be of the same length as sample size")

    roll <- seq(0,1,length.out=n)

    # Separate high and low shuffles into separate piles

    if(n%%2 == 1) warning("rollBlock always allocates a leftover patient to control")

    highRoll <- roll[roll>0.5]
    lowRoll <- roll[roll<=0.5]

    # Shuffle each pile separately
    highRoll <- highRoll[order(runif(length(highRoll)))]
    lowRoll <- lowRoll[order(runif(length(lowRoll)))]

    # Riffle shuffle these together to distribute high and low rolls throughout
    roll <- rep(NA, n)

    roll[2*(1:length(lowRoll))-1] <- lowRoll
    roll[2*(1:length(highRoll))] <- highRoll

    # Take each roll block and shuffle it

    for(i in unique(rollBlocks))
    {
     roll[rollBlocks == i] <- roll[rollBlocks == i][order(runif(sum(rollBlocks == i)))]
    }
  }



  nStrata <- length(unique(stratum))

  rand_out <- csmsb:::pocock_randomisation(data_cleaned,
                                           stratum = stratum_cleaned, nStrata = nStrata,
                                           roll,nClasses = nClasses,biasStrength = biasStrength,
                                           balanceGroups = balanceGroups)
  arm <- factor(rand_out$arm,levels=0:1,labels=c("control","treatment"))

  out <- list(sim_df=data,
              stratum=stratum,
              arm=arm,
              roll=roll,
              prob=rand_out$prob,
              direction=rand_out$direction
  )


  # Look at each interrim plus end of trial
  interrims <- c(interrims,n)

  # Get overall and per-stratum balance at each look
  stratum_balance <- look <- allocationBalance <- expected_bias <- prop_random <- imbalance_std <- imbalance_genodds <- {}

  for(i in interrims){

    stratum_balance <- c(stratum_balance,"Overall")
    look <- c(look, i)

    imbalance_std <- rbind(imbalance_std,getImbalance(data[1:i,],arm[1:i],standard=T))
    imbalance_genodds <- rbind(imbalance_genodds,getImbalance(data[1:i,],arm[1:i],standard = F))
    prop_random <- c(prop_random,getProportionRandom(rand_out$prob[1:i]))
    expected_bias <- c(expected_bias,getExpectedBias(rand_out$prob[1:i]))
    allocationBalance <- rbind(allocationBalance,t(as.matrix(table(arm[1:i]))))


    for(j in unique(stratum)){
      stratum_balance <- c(stratum_balance,j)
      look <- c(look, i)

      imbalance_std <- rbind(imbalance_std,getImbalance(data[1:i,][stratum[1:i]==j,],arm[1:i][stratum[1:i]==j],standard=T))
      imbalance_genodds <- rbind(imbalance_genodds,getImbalance(data[1:i,][stratum[1:i]==j,],arm[1:i][stratum[1:i]==j],standard = F))
      prop_random <- c(prop_random,getProportionRandom(rand_out$prob[1:i][stratum[1:i]==j]))
      expected_bias <- c(expected_bias,getExpectedBias(rand_out$prob[1:i][stratum[1:i]==j]))
      allocationBalance <- rbind(allocationBalance,t(as.matrix(table(arm[1:i][stratum[1:i]==j]))))
    }

  }

  out$look <- look
  out$stratum_balance <- stratum_balance
  out$imbalance_std <- imbalance_std
  out$imbalance_genodds <- imbalance_genodds
  out$prop_random <- prop_random
  out$expected_bias <- expected_bias
  out$allocationBalance <- allocationBalance

  return(out)
}

