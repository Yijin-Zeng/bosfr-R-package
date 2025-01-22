#source('SortingAlg.R')
#source('CountNumbers.R')
#source('UpperBoundCaseOne.R')

# Compute upper bounds of Spearman's footrule under missing case I
ComputeUpperBoundCaseOne <- function(X, Y, fullres = FALSE) {
  # X: potentially partially observed data.
  # Y: fully observed.

  # If fullres = TRUE, returns all (m1 + 1) ranks

  n <- length(X)

  # Reorder X and Y such that missing values in X come first and the correspoding samples in Y increase in order
  Res <- SortingAlg(X,Y)
  X <- Res$SortedX
  Y <- Res$SortedY

  m1 <- sum(is.na(X)) # Number of missing values
  U <- 1:m1 # Indexes of missing values

  # Calculate ranks
  ranks_X <- rank(X, na.last = 'keep')
  ranks_Y <- rank(Y)

  # Calculate d
  d <- ranks_X - ranks_Y
  d <- d[!is.na(d)]

  # Count numbers in d between -m1 to - 1 (as in r), and the total numbers in d larger or equal than 0 (as in R_Positive)
  Res <- CountNumbers(d, seq(-m1,-1))
  r <- Res$Counted_Numbers
  R_Positive <- sum(d >= 0)

  # Calculate initial Spearman footrule distance
  Initial_D <- sum(abs(seq(n, (n-m1+1), -1) - ranks_Y[1:m1] )) + sum(abs(d))

  # Iteratively calculate upper bound for each t
  # for (t in 1:m1) {
  #   D[t+1] <- D[t] + abs(t - ranks_Y[(m1 - t + 1)]) - abs((n - m1 + t) - ranks_Y[(m1 - t + 1)]) + R_Positive - ((n - m1) - R_Positive)
  #   R_Positive <- R_Positive + as.numeric(r[as.character(-t)])
  # }
  Res <- UpperBoundCaseOne(Initial_D = Initial_D, ranks_Y_U = ranks_Y[U], R_Positive = R_Positive, r = r, n = n, fullres = fullres)

  # Return the maximum possible Spearman footrule distance
  return(Res)

}

# ### Example usage:
# source('ElaborateSpearmanFootrule.R')  ## The function to elaborate all possible Speamrna footrule
#
# # Example One
# X <- c(1, 2, 3, NA, 5)
# Y <- c(1, 4, 2, 6, 3)
# ComputeUpperBoundCaseOne(X, Y)
# ElaborateSpearmanFootrule(X,Y)
#
# # Example Two
# X <- c(NA, 2, 3, NA, 5)
# Y <- c(1, 4, 2, 6, 3)
# ComputeUpperBoundCaseOne(X, Y)
# ElaborateSpearmanFootrule(X,Y)
#
# # Example Three
# X <- c(NA, 2, 3, NA, 5, NA, NA,1,10,7)
# Y <- c(1, 4, 2, 6, 3,10,7,8, 9,15)
# ComputeUpperBoundCaseOne(X, Y)
# res1 <- ElaborateSpearmanFootrule(X,Y)
# max(res1$AllPossibleSPearmanFootrule)

