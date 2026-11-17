pocock <- function(D = c("var","sd","range"),
                   D_upper_limit = 0,
                   D_upper_limit_value = c("1","D"),
                   G_weights = 1,
                   p = c("first","rank","value"),
                   p_constant = 0.5
                ){


}

## CPP Implementations ########################################################


pocock_allocate_cpp <- function(data,
                               arm,
                               newData,
                               settings,
                               randSequence,
                               offset
){



}

csmsb_trialSim_cpp <- function(data,
                               settings,
                               randSequence,
                               replay){
    
}

## R Implementations ########################################################
pocock_allocate_R <- function(data,
                              arm,
                              newData,
                              settings,
                              randSequence,
                              offset
){

  # Extract out and process inputs
  rand_number <- nrow(data)+1+offset

  D <- settings$D
  D_upper_limit = settings$D_upper_limit
  D_upper_limit_value = as.character(settings$D_upper_limit_value)
  G_weights = settings$G_weights
  p = setting$p
  p_constant = settings$p_constant

  if(length(G_weights)==1) G_weights <- rep(G_weights,ncol(data))  

  if(!(D_upper_limit_value %in% c("1","D"))) stop("D_upper_limit should be either '1' or 'D'")

  #   pocock and simon assumes everything is a factor
  if(!is.factor(arm)) arm <- as.factor(arm)

  n_groups <- levels(arm)
  
  for(i in colnames(data)){
      data[,i] <- as.factor(data[,i])
      newData[,i] <- factor(newData[,i], levels=levels(data[,i]))
  }


  # Reference paper:  
  #  https://www.jstor.org/stable/2529712  
    
  # Get resultant lack of balance in variable i if patient is assigned to arm k  
  D_ik <- do.call("rbind",lapply(1:ncol(data),function(i){
      sapply(levels(arm),function(k){
          
          # This, if treatment k were assinged to the next patient
          # $d_ik = D({x_{i r_i l}^k}_{l=1}^k)$
          # would be the resultant lack of balance among treatment
          # assignments for patients with level $r_i$ of factor $i$
        
          xtab <- table(c(arm[data[,i] == newData[,i]],k))
          xtab <- c(xtab)

          out <- switch(D,
                        var = var(xtab),
                        sd = sd(xtab),
                        range = diff(range(xtab))
          )

          if(D_upper_limit > 0){
              
            this_D <- out 
            if(D_upper_limit_value == "1") this_D <- 1
            
            if(out<=D_upper_limit){
                out <- 0
            } else {
                out <- this_D
            }
          }
      })
  }))

  G_k <- sapply(1:ncol(D_ik), function(k){
    sum(G_weights * D[i,])  
  })
  
  G_k_order <- order(G_k)

  prob <- rep(NA,length(G_k))

  if(p=="first"){

    prob[G_k_order == 1] <- p_constant
    prob[G_k_order != 1] <- (1-p_constant)/(length(G_k)-1)

  } else if(p=="rank"){

    prob[G_k_order] <- sapply(1:length(G_k),function(k){p_constant - 2*(length(G_k)*p_constant - 1)/(length(G-k)*(length(G_k)+1))*k})
      
  } else if(p=="value"){

    prob <- 1/(length(G_k) - p_constant) * (1 - p_constant * G_k / sum(G_k))

  } 
  
  details <- list(
    D_ik,
    G_k,
    prob = prob
  )


  ## Random Assignment ##################################################

  # Get the random number for this allocation

  if("Runif" %in% class(randSequence$sequence)){

    randSequence <- randSequence$sequence
    roll <- randSequence[rand_number]

    # 0 indexed groups
    decision <- min(which(cumsum(prob)<roll)) - 1

  } else if ("Permutedblock" %in% class(randSequence$sequence)){

    if(p != "first") stop("Permutedblock is only supported when p is \"first\"")

    randSequence <- randSequence$sequence
    roll <- randSequence[rand_number,]

    if(roll[2]==0){
      # non-flexible
      decision <- roll[1] 

    } else {
      # Flexible
                 
      decision <- which(prob==max(prob))-1
      if(length(decision)>1){
        # We have multiple options to choose from.

        # If our original decision is an option, choose it.
        if(roll[1] %in% decision){
          decision <- roll[1]
        } else {
          # Else, fall back on simple randomisation
          choice <- min(which( (1:length(decision))/length(decision) < roll[3]))
          decision <- decision[choice]
        }

      } 

    }

    roll <- paste(roll, collapse = ", ")

  }


  return(
    list(
    decision = decision,
    prob = prob,
    roll = roll,
    biasDirection = biasDirection,
    details = details
  ))


    
    

}