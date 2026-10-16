test_that("genodds implementation works", {

  x1 <- genodds::genodds(mtcars$mpg,mtcars$am)
  x2 <- allocatR::test_genodds(mtcars$mpg,mtcars$am)

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

    df <- mtcars[,c("mpg","am","gear")]
    df$gear <- factor(df$gear)
    rownames(df) <- NULL
    randSequence = get_randNumbers(nrow(df), trng = s)

    tmp <- allocatR::trialSim(csmsb(alpha=0.5, biasParameters=0.7),
                              data=df,replay = T,
                              randSequence = randSequence,
                              nThreads = 1
                              )

    sapply(1:(nrow(df)-1),function(i){

      allocatR::allocate(csmsb(alpha=0.5, biasParameters=0.7),
                         arm = c(control=0,treatment=1)[tmp$arm[1:i]],
                         data=df[1:i,],
                         newData =df[i+1,],
                         randSequence = randSequence,
                         implementation = "R" )$prob

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

  s <- 1

  randSequence = get_randNumbers(nrow(df)+1, trng = s)

  sapply(seq(0,1,0.01), function(alpha){


    allocatR::allocate(csmsb(alpha=alpha, biasParameters=0.7),
                       arm = arm,
                       data=df,
                       newData =df[1,],
                       randSequence = randSequence,
                       implementation = "R")$prob

  }) -> R_implementation

  sapply(seq(0,1,0.01), function(alpha){

    allocatR::allocate(csmsb(alpha=alpha, biasParameters=0.7),
                       arm = arm,
                       data=df,
                       newData =df[1,],
                       randSequence = randSequence,
                       implementation = "cpp")$prob

  }) -> cpp_implementation

  expect_true(all(R_implementation==cpp_implementation))

})


test_that("C++ implementation of cs-msb matches R implementation when weights are uneven", {

  set.seed(514327)
  seedList <- round(runif(100,-.Machine$integer.max+1,.Machine$integer.max-1))

  sapply(seedList,function(s){

    df <- mtcars[,c("mpg","am","gear")]
    df$gear <- factor(df$gear)
    rownames(df) <- NULL
    randSequence = get_randNumbers(nrow(df), trng = s)

    tmp <- allocatR::trialSim( csmsb(alpha=0.5,
                                              biasParameters=0.7,
                                              weights = c(2,1,1)
                                              ),
                              data=df,replay = T,
                              randSequence = randSequence,
                              nThreads = 1
    )

    sapply(1:(nrow(df)-1),function(i){

      allocatR::allocate(csmsb(alpha=0.5, biasParameters=0.7,weights = c(2,1,1)),
                         arm = c(control=0,treatment=1)[tmp$arm[1:i]],
                         data=df[1:i,],
                         newData =df[i+1,],
                         randSequence = randSequence,
                         implementation = "R" )$prob

    }) -> tmp2

    cbind(tmp$probabilityLog[-1],tmp2,abs(tmp$probabilityLog[-1]!=tmp2))

    max(abs(tmp$probabilityLog[-1]!=tmp2))

  }) -> difference_seen

  if(max(difference_seen)>0)
  {
    warning(sprintf("Errors with seeds %s",paste(seedList[difference_seen!=0], collapse=", ")))
  }

  expect_equal(max(difference_seen),0)

})
