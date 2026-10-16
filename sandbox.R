















rm(list=ls())

library(allocatR)

s <- -1914854074

df <- mtcars[,c("mpg","am","gear")]
df$gear <- factor(df$gear)
rownames(df) <- NULL
randSequence = get_randNumbers(nrow(df), trng = s)

tmp <- allocatR::trialSim(msb(alpha=0.5, coinStrength =0.7),
                          data=df,replay = T,
                          randSequence = randSequence,
                          nThreads = 1
                          )

sapply(1:(nrow(df)-1),function(i){

  allocatR::allocate(msb(alpha=0.5, coinStrength =0.7),
                      arm = c(control=0,treatment=1)[tmp$arm[1:i]],
                      data=df[1:i,],
                      newData =df[i+1,],
                      randSequence = randSequence,
                      implementation = "R" )$prob

}) -> tmp2


which(tmp$probabilityLog[-1] != tmp2)
i <- 28

# These disagree with one another. Issue is with the t-test p-value.

  allocatR::allocate(msb(alpha=0.5, coinStrength =0.7),
                      arm = c(control=0,treatment=1)[tmp$arm[1:i]],
                      data=df[1:i,],
                      newData =df[i+1,],
                      randSequence = randSequence,
                      implementation = "cpp" )

  allocatR::allocate(msb(alpha=0.5, coinStrength =0.7),
                      arm = c(control=0,treatment=1)[tmp$arm[1:i]],
                      data=df[1:i,],
                      newData =df[i+1,],
                      randSequence = randSequence,
                      implementation = "R" )














df <- mtcars[,c("mpg","am","gear")]
arm <- rep(0:1, length.out = 11)
df$gear <- factor(df$gear)

randSequence = get_randNumbers(100,n.trials =2 , trng = 1)

allocate(msb(alpha=0.5, coinStrength =0.7),
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         implementation = "R"
)

allocate(msb(alpha=0.5, coinStrength =0.7),
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         implementation = "cpp"
)

allocate(csmsb(alpha=0.5, biasParameters =0.7),
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         implementation = "R"
)




trialSim(csmsb(alpha=0.5, biasParameters =0.7),
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df
)








randSequence = get_permutedBlock(100,blockSize = c(4,6),n.trials =1 , trng = 1)


msb(alpha=0.3,coinStrength = 0.7)


allocate(msb(alpha=0.5, coinStrength =0.7),
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         implementation = "cpp"
)


allocate(csmsb(alpha=0.5, biasParameters=0.7),
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         implementation = "cpp"
)

trialSim(csmsb(alpha=0.5, biasParameters=0.7),
         data = df,
         randSequence = randSequence,
         nThreads = 1
)


allocate(csmsb(alpha=0.5, biasParameters=0.7),
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         implementation = "R"
)




