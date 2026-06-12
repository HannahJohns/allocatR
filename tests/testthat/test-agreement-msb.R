
test_that("C++ implementation of msb matches R implementation", {

  set.seed(514327)
  seedList <- round(runif(100,-.Machine$integer.max+1,.Machine$integer.max-1))
  sapply(seedList,function(s){

    set.seed(s)
    df <- mtcars[,c("mpg","am","gear")]
    df$gear <- factor(df$gear)
    rownames(df) <- NULL
    roll <- runif(nrow(df))
    tmp <- csmsb::trialSim(data=df,replay = T,
                           alpha = 0.9999,method = "MSB",
                           biasParameters = 0.6)


    # cbind(tmp$randProb,tmp$testStatistic)

    # Need to test each randomisation decision in the simulation separately to make sure it agrees with rand decision
    # (i.e. no errors in the calculation of concordant pairs)

    sapply(1:(nrow(df)-1),function(i){
      csmsb::R_decideBias(data=df[1:i,],
                          arm=c(control=0,treatment=1)[tmp$arm[1:i]],
                          newRecord = df[i+1,],
                          method = "MSB",
                          alpha = 0.9999,
                          biasParameters = 0.6
      )$prob
    }) -> tmp2

    tmp$probabilityLog[-1]
    tmp2

    max(abs(tmp$probabilityLog[-1]!=tmp2))
  }) -> difference_seen

  if(max(difference_seen)>0)
  {
    warning(sprintf("Errors with seeds %s",paste(seedList[difference_seen!=0], collapse=", ")))
  }

  expect_equal(max(difference_seen),0)

})
