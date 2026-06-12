test_that("genodds implementation works", {

  x1 <- genodds::genodds(mtcars$mpg,mtcars$am)
  x2 <- test_genodds(mtcars$mpg,mtcars$am)

  apply(rbind(
    c(log(x2[1]),x1$pooled_lnodds),
    c(x2[2]/x2[1],x1$pooled_SElnodds),
    c(x2[3],x1$pooled_SElnnull)
  ),1,diff) -> difference

  expect_true(prod(abs(difference) < 1e-8)==1)
})

test_that("C++ implementation of cs-msb matches R implementation", {

  set.seed(514327)
  seedList <- round(runif(100,-.Machine$integer.max+1,.Machine$integer.max-1))

  sapply(seedList,function(s){

    set.seed(s)
    df <- mtcars[,c("mpg","am","gear")]
    df$gear <- factor(df$gear)
    rownames(df) <- NULL
    roll <- runif(nrow(df))
    tmp <- csmsb::trialSim(data=df,replay = T,
                           alpha = 0.9999,
                           biasParameters = 0.6)

    # cbind(tmp$randProb,tmp$testStatistic)


    # Need to test each randomisation decision in the simulation separately to make sure it agrees with rand decision
    # (i.e. no errors in the calculation of concordant pairs)

    sapply(1:(nrow(df)-1),function(i){
      csmsb::R_decideBias(data=df[1:i,],
                      arm=c(control=0,treatment=1)[tmp$arm[1:i]],
                      newRecord = df[i+1,],
                      method = "CSMSB",
                      voteMethod = "pooled",
                      biasMethod = "static",
                      alpha = 0.9999,
                      biasParameters = 0.6
                      )$prob
    }) -> tmp2

    max(abs(tmp$probabilityLog[-1]!=tmp2))
  }) -> difference_seen

  if(max(difference_seen)>0)
  {
    warning(sprintf("Errors with seeds %s",paste(seedList[difference_seen!=0], collapse=", ")))
  }

  expect_equal(max(difference_seen),0)

})



test_that("C++ and R implementations agree on when threshold for covariates is reached", {


  df <- mtcars[,c("mpg","am","gear")]
  df$gear <- factor(df$gear)
  rownames(df) <- NULL

  arm <- rep(0:1, length.out=nrow(df))

  sapply(seq(0,1,0.01), function(alpha){
    csmsb::R_decideBias(data=df,
                        arm=arm,
                        newRecord = df[1,],
                        method = "CSMSB",
                        voteMethod = "pooled",
                        biasMethod = "static",
                        alpha = alpha,
                        biasParameters = 0.6
    )$prob
  }) -> R_implementation

  sapply(seq(0,1,0.01), function(alpha){

    csmsb::decideBias(data=df,
                      arm=arm,
                      newRecord = df[1,],
                      method = "CSMSB",
                      voteMethod = "pooled",
                      biasMethod = "static",
                      alpha = alpha,
                      biasParameters = 0.6
    )$prob
  }) -> cpp_implementation


  expect_equal(max(R_implementation-cpp_implementation),0)

})



test_that("C++ implementation of cs-msb matches R implementation when weights are uneven", {

  set.seed(514327)
  seedList <- round(runif(100,-.Machine$integer.max+1,.Machine$integer.max-1))
  sapply(seedList,function(s){

    set.seed(s)
    df <- mtcars[,c("mpg","am","gear")]
    df$gear <- factor(df$gear)
    rownames(df) <- NULL
    roll <- runif(nrow(df))
    tmp <- csmsb::trialSim(data=df,replay = T,
                           weights = c(2,1,1),
                           alpha = 0.9999,
                           biasParameters = 0.6)

    # cbind(tmp$randProb,tmp$testStatistic)


    # Need to test each randomisation decision in the simulation separately to make sure it agrees with rand decision
    # (i.e. no errors in the calculation of concordant pairs)

    sapply(1:(nrow(df)-1),function(i){
      csmsb::R_decideBias(data=df[1:i,],
                          arm=c(control=0,treatment=1)[tmp$arm[1:i]],
                          newRecord = df[i+1,],
                          method = "CSMSB",
                          voteMethod = "pooled",
                          biasMethod = "static",
                          weights = c(2,1,1),
                          alpha = 0.9999,
                          biasParameters = 0.6
      )$prob
    }) -> tmp2

    max(abs(tmp$probabilityLog[-1]!=tmp2))
  }) -> difference_seen

  if(max(difference_seen)>0)
  {
    warning(sprintf("Errors with seeds %s",paste(seedList[difference_seen!=0], collapse=", ")))
  }

  expect_equal(max(difference_seen),0)

})


