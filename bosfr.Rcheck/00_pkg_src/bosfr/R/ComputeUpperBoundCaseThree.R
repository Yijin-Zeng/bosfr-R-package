#source('SortingAlg.R')
#source('CountNumbers.R')
#source('UpperBoundCaseThree.R')

# Compute upper bounds of Spearman's footrule under missing case III
ComputeUpperBoundCaseThree <- function(X, Y, fullres = FALSE) {
  # X, Y: potentially unobserved data, for each component i \in [n],
  # X[i] and Y[i] are either both being observed, or both missing.

  n <- length(X)
  # Reorder X and Y such that missing values positioned in the last
  Res <- SortingAlg(X,Y)
  X <- Res$SortedX
  Y <- Res$SortedY

  m3 <- sum(is.na(X)) # Number of missing values
  Psi <- which(is.na(X)) # Indexes of missing values

  # Calculate ranks
  ranks_X <- rank(X, na.last = 'keep')
  ranks_Y <- rank(Y, na.last = 'keep')

  # Initialize imputation
  ranks_X[Psi] <- seq(n, n-m3+1)
  ranks_Y[Psi] <- seq(1, m3)
  ranks_Y[setdiff(seq(1,n), Psi)] <- ranks_Y[setdiff(seq(1,n), Psi)] + m3

  # Calculate d
  d <- ranks_X - ranks_Y
  d_Both <- d[setdiff(seq(1,n), Psi)] # Where samples in both X and Y are observed
  d_Psi <- d[Psi]

  ### The number of non-negative samples of d within [n]\(U \cup \Phi)
  R_Positive <- sum(d_Both >= 0)
  # Count numbers in d within [n]\Psi between -2*m3 to -1 (as in r)
  Res <- CountNumbers(d_Both, seq( (-2*m3 - 1), -1))
  r <- Res$Counted_Numbers

  # Initialize D
  Initial_D <- sum(abs(d))

  Res <- UpperBoundCaseThree(Initial_D = Initial_D, R_Positive = R_Positive, r = r, m3 = m3, n = n, fullres = fullres)

  return(Res)

}

# ### Example usage:
# source('ElaborateSpearmanFootrule.R')  ## The function to elaborate all possible Speamrna footrule
#
# # Example One
# X <- c(NA, 1, NA, 3, 2)
# Y <- c(NA, 4, NA, 2, 1)
# ComputeUpperBoundCaseThree(X, Y)
# Res <- ElaborateSpearmanFootrule(X,Y)
# max(Res$AllPossibleSPearmanFootrule)
#
#
# # Example Two
# X <- c(NA, 1, NA, 3, 2, 7, 9, NA)
# Y <- c(NA, 4, NA, 2, 6, 3, 1, NA)
# ComputeUpperBoundCaseThree(X, Y)
# Res <- ElaborateSpearmanFootrule(X,Y)  # This may take a while
# max(Res$AllPossibleSPearmanFootrule)
#
#
# # Example Three
# X <- c(4,1,9,2,3,NA)
# Y <- c(1,2,5,7,9,NA)
# ComputeUpperBoundCaseThree(X, Y)
# Res <- ElaborateSpearmanFootrule(X,Y)
# max(Res$AllPossibleSPearmanFootrule)
#
#
# # Example Four
# X <- c(4,NA,NA,NA,2,1,9,8)
# Y <- c(1,NA,NA,NA,7,6,9,5)
# ComputeUpperBoundCaseThree(X, Y)
# Res <- ElaborateSpearmanFootrule(X,Y)   # This may take a while
# max(Res$AllPossibleSPearmanFootrule)


