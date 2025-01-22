# Compute lower bounds of Spearman's footrule under missing case I
ComputeLowerBoundCaseOne <- function(X, Y) {
  # X: potentially incomplete data.
  # Y: fully observed vector.

  n <- length(X)
  U <- which(is.na(X)) # indexes of missing values in X

  # Calculate ranks
  ranks_X <- rank(X, na.last = 'keep')
  ranks_Y <- rank(Y)

  while(length(U) > 0) {
    # Find the index in U for which Y has the smallest rank
    minRankIndexInU <- U[which.min(ranks_Y[U])]

    # Assign a real number to X[minRankIndexInU] such that it has the same rank as Y[minRankIndexInU]. Notice that we can assign rank directly
    ranks_X[minRankIndexInU] <- ranks_Y[minRankIndexInU]

    # Update ranks in [n] \ U
    for (i in setdiff(seq(1,n), U)){
      if( ranks_X[i] >=  ranks_Y[minRankIndexInU] ){
        ranks_X[i] <- ranks_X[i] + 1
      }
    }

    # Update U by removing the index just used
    U <- setdiff(U, minRankIndexInU)
  }

  # Calculate the Spearman's footrule distance
  footrule_distance <- sum(abs(ranks_X - ranks_Y))

  return(list(Minfootrule = footrule_distance, Ranks_ImputedX = ranks_X, Ranks_ImputedY = ranks_Y))
}

# Example usage
#X <- c(NA, NA, 3, 4, 5) # X vector with missing values as NA
#Y <- c(3, 2, 1, 4, 5) # Fully observed Y vector
#res <- ComputeLowerBoundCaseOne(X, Y)
#res$Minfootrule
#source('ElaborateSpearmanFootrule.R')
#ElaborateSpearmanFootrule(X,Y)

