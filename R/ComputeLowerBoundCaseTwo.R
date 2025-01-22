# Compute lower bounds of Spearman footrule under missing case II

ComputeLowerBoundCaseTwo <- function(X, Y) {
  # X,Y: potentially incomplete data, for each component i \in [n],
  # at least one of X[i] or Y[i] is observed.

  n <- length(X)

  # Initialize U and Phi sets
  U <- which(is.na(X))
  Phi <- which(is.na(Y))

  # Calculate ranks
  ranks_X <- rank(X, na.last = 'keep')
  ranks_Y <- rank(Y, na.last = 'keep')

  # While both U and Phi have elements
  while( (length(U) > 0) && (length(Phi) > 0) ) {
    # Find the element with the minimum rank in Phi and U
    l1 <- Phi[which.min(ranks_X[Phi])]
    l2 <- U[which.min(ranks_Y[U])]

    # Determine which set to update based on rank comparison
    if(ranks_X[l1] <= ranks_Y[l2]) {
      # Impute Y[l1] to match rank with X
      ranks_Y[l1] <- ranks_X[l1] # Impute

      # Update ranks after imputation
      for (i in setdiff(seq(1,n), Phi)){
        if( ranks_Y[i] >=   ranks_X[l1] ){
          ranks_Y[i] <- ranks_Y[i] + 1
        }
      }

      # Update Phi by removing the index just used
      Phi <- setdiff(Phi, l1)
    } else {
      # Impute X[l2] to match rank with Y
      ranks_X[l2] <- ranks_Y[l2] # Impute

      # Update ranks after imputation
      for (i in setdiff(seq(1,n), U)){
        if( ranks_X[i] >= ranks_Y[l2] ){
          ranks_X[i] <- ranks_X[i] + 1
        }
      }

      # Update U by removing the index just used
      U <- setdiff(U, l2)

    }
  }

  # If U or Phi is not empty, run Algorithm 1
  if(length(U) > 0) {
    Res <- ComputeLowerBoundCaseOne(ranks_X, ranks_Y)
  }else if (length(Phi) > 0){
    Res_Y_X <- ComputeLowerBoundCaseOne(ranks_Y, ranks_X)
    Res <- list(Minfootrule = Res_Y_X$Minfootrule, Ranks_ImputedX = Res_Y_X$Ranks_ImputedY, Ranks_ImputedY = Res_Y_X$Ranks_ImputedX)
  }else{
    footrule_distance <- sum(abs(ranks_X - ranks_Y))
    Res <- list(Minfootrule = footrule_distance, Ranks_ImputedX = ranks_X, Ranks_ImputedY = ranks_Y)
  }

  # Return
  return(Res)
}


#### Examples

# source('ElaborateSpearmanFootrule.R')  ## The function to elaborate all possible Speamrna footrule
#
# # Example one
# X <- c(NA, 2, 3)
# Y <- c(1, NA, 2)
# Res_1 <- ComputeLowerBoundCaseTwo(X, Y)
# Res_1
# ElaborateSpearmanFootrule(X,Y)
#
#
# # Example Two
# X <- c(4, NA, 6, NA)
# Y <- c(NA, 5, NA, 7)
# Res_2 <- ComputeLowerBoundCaseTwo(X, Y)
# Res_2
# ElaborateSpearmanFootrule(X,Y)
#
#
# # Example Three
# X <- c(4, 6, 2, 3, 5, NA)
# Y <- c(1, 3, 2, 4, 7, 6)
# Res_3 <- ComputeLowerBoundCaseTwo(X, Y)
# Res_3
# ElaborateSpearmanFootrule(X,Y)
#
#
# # Example Four
# X <- c(NA, 4, 5, NA, 3, NA)
# Y <- c(2, NA, NA, 4, NA, 1)
# Res_4 <- ComputeLowerBoundCaseTwo(X, Y)
# Res_4
# temp <- ElaborateSpearmanFootrule(X,Y)
# min(temp$AllPossibleSPearmanFootrule)

