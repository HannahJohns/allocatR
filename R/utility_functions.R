allowed_engines <- function(){
  c("Rcpp_yarn2",
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
}


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


