test_that("Results of trialOC are unaffected by thread count", {


  df <- mtcars[,c("mpg","am","gear")]
  df$gear <- factor(df$gear)
  rownames(df) <- NULL

  singeThread <- trialOC(df,
          replay=T,
          alpha=0.25,
          method = "CSMSB",
          voteMethod = "pooled",
          biasMethod = "static",
          biasParameters = 0.6,
          nSims = 100,
          nThreads = 1,
          trng = 354)

  multiThread <- trialOC(df,
          replay=T,
          alpha=0.25,
          method = "CSMSB",
          voteMethod = "pooled",
          biasMethod = "static",
          biasParameters = 0.6,
          nSims = 100,
          nThreads = 2,
          trng = 354)

  expect_true(all( singeThread$prop_random==multiThread$prop_random))

  expect_true(all(singeThread$std_imbalance_mpg==multiThread$std_imbalance_mpg))

})




test_that("Increasing the weights reduces imbalance", {

  df <- mtcars[,c("mpg","am","gear")]
  df$gear <- factor(df$gear)
  rownames(df) <- NULL
  roll <- runif(nrow(df))
  tmp <- csmsb::trialOC(data=df,replay = T,
                         weights = c(1,1,1),
                         alpha = 0.9999,
                         biasParameters = 0.6,
                        nSims = 1000)

  tmp2 <- csmsb::trialOC(data=df,replay = T,
                        weights = c(2,1,1),
                        alpha = 0.9999,
                        biasParameters = 0.6,
                        nSims = 1000)


  tmp3 <- csmsb::trialOC(data=df,replay = T,
                         weights = c(3,1,1),
                         alpha = 0.9999,
                         biasParameters = 0.6,
                         nSims = 1000)

  testthat::expect_all_true(c(
    mean(abs(as.numeric(tmp$genodds_imbalance_mpg))) >  mean(abs(as.numeric(tmp2$genodds_imbalance_mpg))),
    mean(abs(as.numeric(tmp2$genodds_imbalance_mpg))) >  mean(abs(as.numeric(tmp3$genodds_imbalance_mpg)))
  ))

})


