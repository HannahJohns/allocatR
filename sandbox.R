
rm(list=ls())

df <- mtcars[,c("mpg","am","gear")]
arm <- rep(0:1, length.out = 10)
df$gear <- factor(df$gear)

randSequence = get_randNumbers(100,n.trials =2 , trng = 1)
randSequence = get_permutedBlock(100,blockSize = c(4,6),n.trials =2 , trng = 1)

trialSim(method = "csmsb",
         data = df,
         randSequence = randSequence,
         settings = csmsb_settings(alpha=0.5, biasParameters=0.7),nThreads = 1
)

allocate("csmsb",
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         settings = csmsb_settings(alpha=0.5, biasParameters=0.7), implementation = "cpp"
         )

allocate("csmsb",
         arm = arm[1:10],
         randSequence = get_randNumbers(100,n.trials =1 , trng = 1),
         data = df[1:10,], newData = df[11,],
         settings = csmsb_settings(alpha=0.5, biasParameters=0.7), implementation = "R"
)


arm = arm[1:10]
randSequence = get_randNumbers(100,n.trials =1 , trng = 1)
data = df[1:10,]
newData = df[11,]
settings = csmsb_settings(alpha=0.5, biasParameters=0.7)



