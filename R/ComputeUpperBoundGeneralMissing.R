#source('SortingAlg.R')
#source('CountNumbers.R')
#source('ElaborateSpearmanFootrule.R')
#source('UpperBoundCaseOne.R')
#source('UpperBoundCaseTwo.R')
#source('ComputeUpperBoundCaseOne.R')
#source('ComputeUpperBoundCaseTwo.R')
#source('ComputeUpperBoundCaseThree.R')

# Compute upper bounds of Spearman's footrule under general missing case
ComputeUpperBoundGeneralMissing <- function(X, Y, fullres = FALSE) {
  # X,Y: partially observed data

  # If fullres = TRUE, return all possible (m1+1)(m2+1)(m3+1) ranks
  n <- length(X)

  # Identify missing values in X and Y
  U <- which(is.na(X) & !is.na(Y))
  Phi <- which(is.na(Y) & !is.na(X))
  Psi <- which(is.na(X) & is.na(Y))

  # Number of missing samples
  m1 <- length(U)
  m2 <- length(Phi)
  m3 <- length(Psi)

  ## No Missing data
  if((m1 + m2 + m3 == 0)){
    return(sum(abs(rank(X)-rank(Y))))
    }

  ## Missing case I
  if( (m3 == 0) & (m1 != 0) & (m2 == 0) ){
    Res <- ComputeUpperBoundCaseOne(X,Y,fullres)
    return(Res)
  }
  if( (m3 == 0) & (m1 == 0) & (m2 != 0) ){
    Res <- ComputeUpperBoundCaseOne(Y,X,fullres)
    return(Res)
  }

  ## Missing case II
  if( (m3 == 0) & (m1 != 0) & (m2 != 0) ){
    Res <- ComputeUpperBoundCaseTwo(X,Y,fullres)
    return(Res)
  }

  ## Missing case III
  if( (m3 != 0) & (m1 == 0) & (m2 == 0) ){
    Res <- ComputeUpperBoundCaseThree(X,Y,fullres)
    return(Res)
  }

  ########## Missing Case where U \cup Phi \neq \emptyset and Psi \neq \emptyset. Without loss of generality, let us assume U \neq \emptyset.
  if( (m1 == 0) & (m2 != 0) ){
    temp <- X
    X <- Y
    Y <- temp
  }

  # Reorder X and Y such that missing values positioned in the last
  Res <- SortingAlg(X,Y)
  X <- Res$SortedX
  Y <- Res$SortedY

  # Identify missing values in X and Y
  U <- which(is.na(X) & !is.na(Y))
  Phi <- which(is.na(Y) & !is.na(X))
  Psi <- which(is.na(X) & is.na(Y))

  # Number of missing samples
  m1 <- length(U)
  m2 <- length(Phi)
  m3 <- length(Psi)

  # Calculate ranks
  ranks_X <- rank(X, na.last = 'keep')
  ranks_Y <- rank(Y, na.last = 'keep')

  # Initialize imputation
  ranks_X[Psi] <- seq(n, n-m3+1)
  ranks_Y[Psi] <- seq(1, m3)
  ranks_Y[setdiff(seq(1,n), union(Psi,Phi))] <- ranks_Y[setdiff(seq(1,n), union(Psi,Phi))] + m3
  ranks_X[U] <- seq(n-m3, n-m3-m1+1)
  ranks_Y[Phi] <- seq(n, n-m2+1)

  # Calculate d
  d <- ranks_X - ranks_Y
  d_U <- d[U]
  d_Phi <- d[Phi]
  d_Psi <- d[Psi]
  d_Both <- d[setdiff(seq(1,n), union(U, union(Psi, Phi)))] # Where samples in both X and Y are observed

  ########### Initilize parameters

  ### The number of non-negative samples of d within U
  K_Positive <- sum(d_U >= 0)
  K_Negative <- sum(d_U <= 0)
  # Count numbers in d within U between -2*m3-1 to m2 (as in k)
  Res <- CountNumbers(d_U, seq((-2*m3-1), m2))
  k <- Res$Counted_Numbers  ### k denotes the counting number from -2*m3-1 to m2
  k_truncated <- k[(2*m3+2+1):(2*m3+2+m2)]   ### k_truncated denotes the counting number from 1 to m2

  ### The number of non-negative samples of d within Phi
  S_Positive <- sum(d_Phi >= 0)
  # Count numbers in d within Phi between -m1-2*m3 to - 1 (as in s)
  Res <- CountNumbers(d_Phi, seq((-m1-2*m3), -1))
  s <- Res$Counted_Numbers
  s_truncated <- s[(m1+2*m3+1-m1):(m1+2*m3)]  # s_truncated denotes the counting number from -m1 to -1,

  ### The number of non-negative samples of d within [n]\(U \cup \Phi \cup \Psi)
  R_Positive <- sum(d_Both >= 0)
  R_Negative <- sum(d_Both <= 0)
  # Count numbers in d within [n]\(U \cup \Phi) between -2*m3-m1 to m2 (as in r)
  Res <- CountNumbers(d_Both, seq((-2*m3-m1), m2))           ###
  r <- Res$Counted_Numbers
  r_truncated <- r[(2*m3+m1+1-m1):(2*m3+m1+1+m2)] # r_truncated denotes the counting number from -m1 to m2,

  # Initialize D
  D <- array(data = NA, dim = c((m2+1), (m1+1), (m3+1)))

  # Calculate initial Spearman footrule distance
  D[1,1,1] <- sum(abs(d))

  #Res <- UpperBoundCaseOne(Initial_D = ConditionD[1], ranks_Y_U = ranks_Y[U], R_Positive = R_Positive, r = r, n = n, m3 = m3, t3 = 0)
  Res <- UpperBoundCaseTwo(Initial_D = D[1,1,1], ranks_Y_U = ranks_Y[U], ranks_X_Phi = ranks_X[Phi], R_Positive = R_Positive, R_Negative = R_Negative, S_Positive = S_Positive, K_Negative = K_Negative,  r = r_truncated, k = k_truncated, s = s_truncated, n = n, m3 = m3, t3 = 0, fullres = TRUE)
  D[,,1] <- Res

  ## Iteratively update imputations in Psi
  #r_k_s_MinusOne <- r[as.character(-1)] + k[as.character(-1)] + s[as.character(-1)]
  r_k_s_MinusOne <- r[2*m3+m1] + k[2*m3+1] + s[2*m3+m1] # r denotes the counting number from -2*m3-m1 to m2, k denotes -2*m3-1 to m2, s denotes -m1-2*m3 to - 1
  # Iteratively calculate upper bound for each t
  for (t in 1:m3) {

    a <- t - (n + 1 - t)      # the differences of ranks in n - t + 1 after altering
    b <- (n - m3 + t) - (m3 - t + 1) # the differences of ranks in n - t + 1 before altering

    R_K_S_Positive <- R_Positive + K_Positive + S_Positive

    D[1,1,t+1] <- D[1,1,t] + abs(a) - abs(b) + 2*R_K_S_Positive - 2*((n - m3) - R_K_S_Positive) + 2*r_k_s_MinusOne

    ### Update parameters for calculating next Condition D
    # R_Positive <- R_Positive + as.numeric(r[as.character(-2*t + 1)]) + as.numeric(r[as.character(-2*t)])
    # K_Positive <- K_Positive + as.numeric(k[as.character(-2*t + 1)]) + as.numeric(k[as.character(-2*t)])
    # S_Positive <- S_Positive + as.numeric(s[as.character(-2*t + 1)]) + as.numeric(s[as.character(-2*t)])

    R_Positive <- R_Positive + r[2*m3+m1+1 -2*t + 1] + r[2*m3+m1+1 -2*t]# r denotes the counting number from (-2*m3-m1) to m2
    K_Positive <- K_Positive + k[2*m3+2 -2*t + 1] + k[2*m3+2 -2*t] # k denotes the counting number from (-2*m3-1) to m2
    S_Positive <- S_Positive + s[2*m3+m1+1 -2*t + 1] + s[2*m3+m1+1 -2*t] # s denotes the counting number from (-m1-2*m3) to -1

    # r_MinusOne <- r[as.character( (-2*t -1) )]
    # k_MinusOne <- k[as.character( (-2*t -1) )]
    # s_MinusOne <- s[as.character( (-2*t -1) )]

    r_MinusOne <- r[2*m3+m1+1 -2*t - 1]
    k_MinusOne <- k[2*m3+2 -2*t - 1]
    s_MinusOne <- s[2*m3+m1+1 -2*t - 1]

    r_k_s_MinusOne <- r_MinusOne + k_MinusOne + s_MinusOne
    ### Update parameters for calculating MaxD

    # Update parameters in U
    ranks_Y[U] <- ranks_Y[U] - 1
    #K_Negative <- K_Negative - as.numeric(k[as.character(-2*t+1)]) - as.numeric(k[as.character(-2*t+2)])
    K_Negative <- K_Negative - k[2*m3+2 -2*t + 1] - k[2*m3+2 -2*t + 2]

    k_star <- rep(NA, m2)  # k_star denotes the counting number from 1 to m2,
    for (i in seq(1, m2)) {
      # k_star[as.character(i)] <- k[as.character(i - 2*t)]
      k_star[i] <- k[2*m3+2 + i - 2*t]
    }

    # Update parameters in Phi
    ranks_X[Phi] <- ranks_X[Phi] + 1

    s_star <- rep(NA, m1)  # s_star denotes the counting number from -m1 to -1,
    for (i in seq(-m1, -1)) {
      # s_star[as.character(i)] <- s[as.character(i - 2*t)]
      s_star[m1+1+i] <- s[2*m3+m1+1 +i-2*t]
    }

    # Update parameters in [n] \ (U \cup Phi \cup Psi)
    #R_Negative <- R_Negative - as.numeric(r[as.character(-2*t+1)]) - as.numeric(r[as.character(-2*t+2)])
    R_Negative <- R_Negative - r[2*m3+m1+1-2*t+1] - r[2*m3+m1+1-2*t+2] # r denotes the counting number from (-2*m3-m1) to m2

    r_star <- rep(NA, (m2+m1+1))  # r_star denotes the counting number from -m1 to m2,
    for (i in seq(-m1, m2)) {
      #r_star[as.character(i)] <- r[as.character(i - 2*t)]
      r_star[m1+1+i] <- r[2*m3+m1+1+i-2*t]
    }

    Res <- UpperBoundCaseTwo(Initial_D = D[1,1,t+1], ranks_Y_U = ranks_Y[U], ranks_X_Phi = ranks_X[Phi], R_Positive = R_Positive, R_Negative = R_Negative, S_Positive = S_Positive, K_Negative = K_Negative, r = r_star, k = k_star, s = s_star, n = n, m3 = m3, t3 = t, fullres =TRUE)

    D[,,t+1] <- Res

  }


  # Return the maximum possible Spearman footrule distance
  if(fullres == TRUE){
    return(D)
  }else{
    return(max(D))
  }

}

# ### Example usage:
# source('ElaborateSpearmanFootrule.R')  ## The function to elaborate all possible Speamrna footrule
#
# # Example One
#
# X <- c(1,5, NA, NA, 2,4, NA,NA)
# Y <- c(2,7, 6, 1, 3, 5, NA, NA)
# ComputeUpperBoundGeneralMissing(X, Y)
# Res <- ElaborateSpearmanFootrule(X,Y)
# max(Res$AllPossibleSPearmanFootrule)
#
#
# # Example Two
# X <- c(NA, 1, NA, 3, 7, 9, NA)
# Y <- c(5, 4, NA, 2, NA, 6, 8)
# ComputeUpperBoundGeneralMissing(X, Y)
# Res <- ElaborateSpearmanFootrule(X,Y)
# max(Res$AllPossibleSPearmanFootrule)
#
#
# # Example Three
# X <- c(NA, 1, NA, 3, 7, 9, 2)
# Y <- c(5, 4, NA, 2, NA, 6, NA)
# ComputeUpperBoundGeneralMissing(X, Y)
# Res <- ElaborateSpearmanFootrule(X,Y)
# max(Res$AllPossibleSPearmanFootrule)


#X1 <- c(1,2,3,4,5,6,7,8,NA,NA)
#X2 <- c(3,7,6,1,5,2,4,8,NA,NA)
#X3 <- c(8,7,6,5,4,3,2,1,NA,NA)
#X4 <- c(1,2,3,4,5,6,NA,NA,7,8)
#X5 <- c(1,2,3,4,NA,NA,5,6,7,8)
#X6 <- c(1,2,3,4,5,6,7,NA,NA,NA)
#Y <- seq(1,10)
#ComputeUpperBoundGeneralMissing(X1, Y, fullres=TRUE)
#ComputeUpperBoundGeneralMissing(X2, Y, fullres=TRUE)
#ComputeUpperBoundGeneralMissing(X3, Y, fullres=TRUE)
#ComputeUpperBoundGeneralMissing(X4, Y, fullres=TRUE)
#ComputeUpperBoundGeneralMissing(X5, Y, fullres=TRUE)
#ComputeUpperBoundGeneralMissing(X6, Y, fullres=TRUE)
