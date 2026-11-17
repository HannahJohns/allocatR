
K <- 3

p_constant <- (1/K+0.1)

permutations <- gtools::permutations(K,K)

sapply(1:nrow(permutations),function(i){

  Gk_order <- permutations[i,]
  prob <- rep(NA,K)

  prob[Gk_order == 1] <- p_constant
  prob[Gk_order != 1] <- (1-p_constant)/(K-1)

  # prob[Gk_order] <- sapply(1:K,function(k){p_constant - 2*(K*p_constant - 1)/(K*(K+1))*k})

  prob
}) |> t() |> as.data.frame() -> verteces

unique(verteces)

# 1, 2, 3, flexible
c(rep((1-p_constant)/(K-1),K),p_constant-(1-p_constant)/(K-1))


sapply(1:nrow(permutations),function(i){

  Gk_order <- permutations[i,]
  prob <- rep(NA,K)

  # prob[Gk_order == 1] <- p_constant
  # prob[Gk_order != 1] <- (1-p_constant)/(K-1)

  prob[Gk_order] <- sapply(1:K,function(k){p_constant - 2*(K*p_constant - 1)/(K*(K+1))*k})

  prob
}) |> t() |> as.data.frame() -> verteces

unique(verteces)











t(apply(unique(verteces),1,cumsum))


library(tidyverse)
library(prefviz)
tern22 <- as_ternable(verteces)
input_df <- get_tern_data2d(tern22)


p <- ggplot(input_df, aes(x = x1, y = x2)) +
  # Draw the ternary space as an equilateral triangle
  add_ternary_base() + 
  # Plot the observations as points
  geom_point() + 
  # Add vertex labels, taken from the ternable object
  add_vertex_labels(tern22$simplex_vertices)+
  geom_ternary_region(
    x1 = 1/3, x2 = 1/3, x3 = 1/3, # Default reference points. Must sum to 1
    vertex_labels = tern22$vertex_labels, # Labels for the regions
    aes(fill = after_stat(vertex_labels)), 
    alpha = 0.3, color = NA, show.legend = FALSE
  ) 

p









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




