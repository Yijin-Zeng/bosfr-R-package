#source('SortingAlg.R')
#source('UpperBoundCaseOne.R')
#source('CountNumbers.R')
#source('UpperBoundCaseTwo.R')

# Compute upper bounds of Spearman's footrule under missing case II
ComputeUpperBoundCaseTwo <- function(X, Y, fullres = FALSE) {
  # X,Y: potentially incomplete data, for each component i \in [n],
  # at least one of X[i] or Y[i] is observed.

  n <- length(X)

  # Reorder X and Y
  Res <- SortingAlg(X,Y)
  X <- Res$SortedX
  Y <- Res$SortedY

  U <- which(is.na(X)) # Indexes of missing values in X
  m1 <- length(U) # Number of missing values in X
  Phi <- which(is.na(Y)) # Indexes of missing values in Y
  m2 <- length(Phi) # Number of missing values in Y

  # Calculate ranks
  ranks_X <- rank(X, na.last = 'keep')
  ranks_Y <- rank(Y, na.last = 'keep')

  # Initialize imputation
  ranks_X[U] <- seq(n, n-m1+1)
  ranks_Y[Phi] <- seq(n, n-m2+1)

  # Calculate d
  d <- ranks_X - ranks_Y
  d_U <- d[U]
  d_Phi <- d[Phi]
  d_Both <- d[setdiff(seq(1,n), union(U, Phi))] # Where samples in both X and Y are observed

  ### The number of non-negative samples of d within U
  K_Negative <- sum(d_U <= 0)
  # Count numbers in d within U between 0 to m2 - 1 (as in k)
  Res <- CountNumbers(d_U, seq(1, m2))
  k <- Res$Counted_Numbers

  ### The number of non-negative samples of d within Phi
  S_Positive <- sum(d_Phi >= 0)
  # Count numbers in d within Phi between -m1 to - 1 (as in s)
  Res <- CountNumbers(d_Phi, seq(-m1, -1))
  s <- Res$Counted_Numbers

  ### The number of non-negative samples of d within [n]\(U \cup \Phi)
  R_Positive <- sum(d_Both >= 0)
  R_Negative <- sum(d_Both <= 0)
  # Count numbers in d within [n]\(U \cup \Phi) between -m1 to m2-1 (as in r)
  Res <- CountNumbers(d_Both, seq(-m1, m2))
  r <- Res$Counted_Numbers

  Initial_D <- sum(abs(d))

  Res <- UpperBoundCaseTwo(Initial_D = Initial_D, ranks_Y_U = ranks_Y[U], ranks_X_Phi = ranks_X[Phi], R_Positive = R_Positive, R_Negative = R_Negative, S_Positive = S_Positive, K_Negative = K_Negative,  r = r, k = k, s = s, n = n, fullres = fullres)

  # Return the maximum possible Spearman's footrule distance
  return(Res)

}


#
# # ### Example usage:
# source('ElaborateSpearmanFootrule.R')
#
# ## Example One
# X <- c(NA, 1, 5, 2, 3)
# Y <- c(6, NA, NA, 4, 2)
# ComputeUpperBoundCaseTwo(X,Y)
# res1 <- ElaborateSpearmanFootrule(X,Y)
# max(res1$AllPossibleSPearmanFootrule)
#
#
# ## Example Two
# X <- c(NA, 1, 5, NA, NA)
# Y <- c(6, NA, NA, 4, 2)
# ComputeUpperBoundCaseTwo(X,Y)
# res1 <- ElaborateSpearmanFootrule(X,Y)
# max(res1$AllPossibleSPearmanFootrule)
#
#
# ## Example Three
# X <- c(NA, 1)
# Y <- c(1,NA)
# ComputeUpperBoundCaseTwo(X,Y)
# res1 <- ElaborateSpearmanFootrule(X,Y)
# max(res1$AllPossibleSPearmanFootrule)
#
#
# ## Example Four
# X <- c(NA, NA, NA, 2,4,6,7,1, 8)
# Y <- c(1,5,7,9,10, NA, NA, NA, 2)
# ComputeUpperBoundCaseTwo(X,Y)
# res1 <- ElaborateSpearmanFootrule(X,Y) # This may take a while
# max(res1$AllPossibleSPearmanFootrule)
#
#
# ## Example Five
# X <- c(NA, NA, NA, 2,4)
# Y <- c(1,5,7,9,10)
# ComputeUpperBoundCaseTwo(X,Y)
# res1 <- ElaborateSpearmanFootrule(X,Y)
# max(res1$AllPossibleSPearmanFootrule)
