get_randNumbers <- function(n,
                n.trials = 1,
                trng = NULL,
                jump_trng = 0
){

  # RNG preamble

  if(!is.null(trng)){
    if(!("rTRNG" %in% rownames(installed.packages()))){
      stop("rTRNG package must be installed to use trng option")
    }

    # Allow passing of a seed, instead of a fully specified engine
    if(is.numeric(trng)) {

      if(as.integer(trng) != trng){
        stop("trng must be a whole number if a seed is being specified")
      }
      trng <- rTRNG::yarn3$new(trng)

    }

    if(!(class(trng) %in% allowed_engines())){

      stop(sprintf("%s is not an allowed rTRNG parallel engine. trng must be an integer, or any of the following types: %s",
                   class(trng),
                   paste(gsub("Rcpp_","",allowed_engines()), collapse=", ")
      ))

    }


    # Jump seed state forward. Use multiple jumps if neccesary.
    while(jump_trng > .Machine$integer.max){
      trng$jump(.Machine$integer.max)
      jump_trng <- jump_trng - .Machine$integer.max
    }

    if(jump_trng > 0){
      trng$jump(jump_trng)
    }

  }


  if(is.null(trng)){
    roll <- runif(n*n.trials)
  } else {
    roll <- rTRNG::runif_trng(n*n.trials,engine=trng)
  }

  if(n.trials>1)  roll <- matrix(roll, nrow = n, ncol = n.trials)

  roll

}


get_permutedBlock <- function(n,
                              blockSize,
                              prop_flexible=0,
                              n.trials=1,
                              max_prop_flexible = 1,
                              max_attempts = 1000,
                              trng = NULL,
                              jump_trng = 0
                              )
{

  if(length(prop_flexible) == 1){
    if(length(prop_flexible) == 1) prop_flexible <- rep(prop_flexible,length(blockSize))
  }

  if(any(prop_flexible > 0)){

    if(any(((blockSize*prop_flexible) %% 2) != 0)) stop(sprintf("NEED AND ERROR MESSAGE HERE"))

  } else {
    if(any((blockSize %% 2) != 0)) stop(sprintf("Block sizes should be even"))

  }

  # RNG preamble

  if(!is.null(trng)){
    if(!("rTRNG" %in% rownames(installed.packages()))){
      stop("rTRNG package must be installed to use trng option")
    }

    # Allow passing of a seed, instead of a fully specified engine
    if(is.numeric(trng)) {

      if(as.integer(trng) != trng){
        stop("trng must be a whole number if a seed is being specified")
      }
      trng <- rTRNG::yarn3$new(trng)

    }

    if(!(class(trng) %in% allowed_engines())){

      stop(sprintf("%s is not an allowed rTRNG parallel engine. trng must be an integer, or any of the following types: %s",
                   class(trng),
                   paste(gsub("Rcpp_","",allowed_engines()), collapse=", ")
      ))

    }


    # Jump seed state forward. Use multiple jumps if neccesary.
    while(jump_trng > .Machine$integer.max){
      trng$jump(.Machine$integer.max)
      jump_trng <- jump_trng - .Machine$integer.max
    }

    if(jump_trng > 0){
      trng$jump(jump_trng)
    }

  }

  # For each trial
  out <- lapply(1:n.trials, function(i_trial){


    make_attempt <- T
    n_attempts <- 0
    while(make_attempt){

      n_attempts <- n_attempts + 1

      # Create a random sequence for block size.
      # Ideally we'd use parts() from partition to do this but practically speaking
      # the sample size is going to be far too big to do this properly

      # Approach taken here is:
      # 1. Overshoot
      # 2. Cut as many blocks as needed to get under sample size
      # 3. Regenerate this using only the smaller blocks
      # 4. Shuffle the complete block size order.

      # This induces some bias in block size towards smaller block sizes, but
      # this shouldn't matter in practical terms.

      # blockSize = c(4,6,8)
      # prop_flexible <- c(0.5,1/3,1/4)
      # n <- 52

      blockSizeSeq <- {}
      allowed_Sizes <- rep(TRUE, length(blockSize))
      while(sum(blockSizeSeq[,1]) < n & sum(allowed_Sizes)>0){

        # Add random blocks, even if it's overshooting the target sample size
        while(sum(blockSizeSeq[,1]) < n){

          if(is.null(trng)){
            roll <- runif(1)
          } else {
            roll <- rTRNG::runif_trng(1,engine=trng)
          }
          roll <- ceiling(roll*sum(allowed_Sizes))

          blockSizeSeq <- rbind(blockSizeSeq, c( blockSize[allowed_Sizes][roll], prop_flexible[allowed_Sizes][roll] ) )

        }

        # Reverse backwards to correct for overshoot.
        while(sum(blockSizeSeq[,1]) > n){
          blockSizeSeq <- blockSizeSeq[1:(nrow(blockSizeSeq)-1),]
        }

        # Eliminate block sizes that are too big
        allowed_Sizes[blockSize > n-sum(blockSizeSeq[,1])] <- FALSE

      }

      if(sum(blockSizeSeq[,1]) == n){

        # Shuffle the blocks around
        if(is.null(trng)){
          roll <- runif(nrow(blockSizeSeq))
        } else {
          roll <- rTRNG::runif_trng(nrow(blockSizeSeq),engine=trng)
        }

        blockSizeSeq <- blockSizeSeq[order(roll),]


      } else {

        warning(sprintf("Failed to generate a block sequence of exactly %d assignments. Final block has been trimmed to size.",n))

        # Exhausted all options and were unable to hit the target number.
        # Take whatever we generated, shuffle, and then append whatver the smallest number is to the end.

        # Shuffle the blocks around
        if(is.null(trng)){
          roll <- runif(nrow(blockSizeSeq))
        } else {
          roll <- rTRNG::runif_trng(nrow(blockSizeSeq),engine=trng)
        }

        blockSizeSeq <- blockSizeSeq[order(roll),]

        allowed_Sizes[blockSize == min(blockSize)] <- TRUE

        if(is.null(trng)){
          roll <- runif(1)
        } else {
          roll <- rTRNG::runif_trng(1,engine=trng)
        }
        roll <- ceiling(roll*sum(allowed_Sizes))

        blockSizeSeq <- rbind(blockSizeSeq, c( blockSize[allowed_Sizes][roll], prop_flexible[allowed_Sizes][roll] ) )

      }

      # For each block sequence, generate the permuted block

      rand_sequence <- do.call("rbind",lapply(1:nrow(blockSizeSeq),function(i){

        thisBlock <- cbind(rep(0:1,length.out = blockSizeSeq[i,1]), rep(0,length.out = blockSizeSeq[i,1]))

        # Number of flexible assignments per arm
        n_flex <- (blockSizeSeq[i,1] * blockSizeSeq[i,2])/2

        if(n_flex >0){
          thisBlock[thisBlock[,1]==0,2][1:n_flex] <- 1
          thisBlock[thisBlock[,1]==1,2][1:n_flex] <- 1
        }

        # Shuffle the blocks around
        if(is.null(trng)){
          roll <- runif(nrow(thisBlock))
        } else {
          roll <- rTRNG::runif_trng(nrow(thisBlock),engine=trng)
        }

        thisBlock <- thisBlock[order(roll),]
        thisBlock

      }))

      rand_sequence <- rand_sequence[1:n,]

      if(mean(rand_sequence[,2]) <= max_prop_flexible){
        make_attempt <- F
      } else {
        if(n_attempts > max_attempts) stop("Was unable to generate a sequence with proportion of flexible assignments below threshold.")
      }
    }


    rand_sequence

  })

  if(n.trials == 1){
    out <- out[[1]]
  } else {
    out
  }

  out

}
