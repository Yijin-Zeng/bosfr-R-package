# Compute lower bounds of Spearman's footrule under missing case III

ComputeLowerBoundCaseThree <- function(X, Y) {
  # X, Y: potentially unobserved data, for each component i \in [n],
  # X[i] and Y[i] are either both being observed, or both missing.
  n <- length(X)
  ranks_X <- rank(X, na.last = NA)
  ranks_Y <- rank(Y, na.last = NA)
  Res <- sum(abs(ranks_X - ranks_Y))
  # Return
  return(list(Minfootrule = Res))
}
