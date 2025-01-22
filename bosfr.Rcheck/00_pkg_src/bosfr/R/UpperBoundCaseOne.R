# Compute the upper bound of Spearman's footrule given:
# (1): Initial_D, the value of Spearman's footrule when no data within U smaller than observed data in X
# (2): ranks_Y_U, where U = {1,2,...,m1}, a vector including the ranks of data in Y within U.
# (3): R_Positive, the number of samples within [n] \ U such that the ranks of X larger or equal than its paired samples
# (4): r, a vector including the counting number in d between -m1 to - 1
# (5): n, sample size

UpperBoundCaseOne <- function(Initial_D, ranks_Y_U, R_Positive, r, n, m3 = 0, t3 = 0, fullres = FALSE) {
  # m3, t3 considers the case where Psi \neq \emptyset, thus it is not missing case I
  # m3 denotes the size of Psi,
  # t3 denotes the number of samples in Psi, smaller than all the other samples in [n] \ Psi

  # If fullres = TRUE, returns all (m1+1) results

  # Number of unobserved samples in X
  m1 <- length(ranks_Y_U)

  D <- rep(0, (m1 + 1))

  # Initial Spearman footrule distance
  D[1] <- Initial_D

  # Iteratively calculate upper bound for each t
  for (t in 1:m1) {
    a <- t3 + t - ranks_Y_U[(m1 - t + 1)]
    b <- (n - m1 + t) - (m3 - t3) - ranks_Y_U[(m1 - t + 1)]
    D[t+1] <- D[t] + abs(a) - abs(b) + R_Positive - (n - m1 - m3 - R_Positive)
    # R_Positive <- R_Positive + as.numeric(r[as.character(-t)])
    R_Positive <- R_Positive + r[-t + m1 + 1] # r denotes the counting number from -m1 to 1,
  }

  # Return the maximum possible Spearman footrule distance
  if(fullres == TRUE){
    return(D)
  }else{
    return(max(D))
  }

}




