
rm(list=ls())

library(allocatR)

# r <- yarn3$new(1)
#
# r1 <- yarn3$new(r$toString())
#
# all(runif_trng(10,engine=r)==runif_trng(10,engine=r1))

data <- synthetic_covariate_data

randScheme <- simple()
randSeq2 <- get_permutedBlock(100, blockSize = 10,n.trials = 2, trng=1)
allocatR::trialSim(randomisationScheme=randScheme, randSequence=randSeq2,  data=data, replay=FALSE, nThreads = 1)



test <- function(randSequence){
  randSeq2$sequence[[1]]
}

test(randSeq2)


randSeq <- get_randNumbers(100,n.trials = 2, trng=1)
trialSim(randomisationScheme=randScheme, randSequence=randSeq,  data=data, replay=FALSE, nThreads = 1)



randScheme <- csmsb(alpha = 0.5, biasParameters=0.7)

trialSim(randomisationScheme=randScheme, randSequence=randSeq,  data=data, replay=FALSE)

randSeq2 <- get_randNumbers(100, n.trials=10, trng=2)

trialSim(randomisationScheme=randScheme, randSequence=randSeq2,  data=data, replay=FALSE)

class(randSeq2$trng)

df <- mtcars[,c("mpg","am","gear")]
arm <- rep(0:1, length.out = 10)
df$gear <- factor(df$gear)

randSequence = get_randNumbers(100,n.trials =2 , trng = 1)
randSequence = get_permutedBlock(100,blockSize = c(4,6),n.trials =1 , trng = 1)


trialSim(csmsb(alpha=0.5, biasParameters=0.7),
         data = df,
         randSequence = randSequence,
         nThreads = 1
)



allocate(csmsb(alpha=0.5, biasParameters=0.7),
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         implementation = "cpp"
         )

allocate(csmsb(alpha=0.5, biasParameters=0.7),
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         implementation = "R"
)



