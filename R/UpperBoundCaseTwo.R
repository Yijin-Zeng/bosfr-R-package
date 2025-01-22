# Compute the upper bound of Spearman's footrule under missing case II given:
# Initial_D, ranks_Y_U, ranks_X_Phi, R_Positive, R_Negative, S_Positive, K_Negative,  r, k, s, n

UpperBoundCaseTwo <- function(Initial_D, ranks_Y_U, ranks_X_Phi, R_Positive, R_Negative, S_Positive, K_Negative,  r, k, s, n, m3 = 0, t3 = 0, fullres = FALSE) {
  # m3, t3 considers the case where Psi \neq \emptyset.
  # m3 denotes the size of Psi,
  # t3 denotes the number of samples in Psi, smaller than all the other samples in [n] \ Psi

  m1 = length(ranks_Y_U)
  m2 = length(ranks_X_Phi)

  # Initialize D
  D <- matrix(data = 0, nrow = (m2 + 1), ncol = (m1+1))


  # Initial Spearman footrule distance
  D[1,1] <- Initial_D

  ####### Calculate the maximum Spearman footrule given initial imputation of Y
  r_s <- rep(0, m1)
  # names(r_s) <- as.character(seq(-m1,-1))
  #
  # for (i in seq(-m1, -1)) {
  #   r_s[as.character(i)] <- r[as.character(i)] + s[as.character(i)]
  # }


  for (i in seq(1,m1)) {
    r_s[i] <- r[i] + s[i]   # r denotes the counting number from -m1 to m2, s denotes the counting number from -m1 to -1,
    # r_s denotes the counting number from -m1 to -1
  }

  R_S_Positive <- R_Positive + S_Positive

  Res <- UpperBoundCaseOne(Initial_D = D[1,1], ranks_Y_U = ranks_Y_U, R_Positive = R_S_Positive, r = r_s, n = n, m3 = m3, t3 = t3, fullres = TRUE)
  D[1,] <- Res

  ##### Iteratively update imputations in Y
  if(m2 == 0){  ## Case I
    if(fullres == TRUE){
     return(D)
    }else{
      return(max(D))
    }
  }

  for(t in 1:m2){

    a <- (m3 - t3) + t - ranks_X_Phi[(m2 - t + 1)]      # the differences of ranks in m1 + m2 - t +1 after altering
    b <- (n - m2 + t) - t3 - ranks_X_Phi[(m2 - t + 1)]  # the differences of ranks in m1 + m2 - t +1 before altering

    D[t+1,1] <- D[t,1] + abs(a) - abs(b) + R_Negative - ((n-m1-m2-m3) - R_Negative) + K_Negative - (m1 - K_Negative)

    ### Update Parameters for updating Spearman footrule of different imputations in Y givn initial imputations in X
    # R_Negative <- R_Negative + r[as.character(t)]
    # K_Negative <- K_Negative + k[as.character(t)]

    R_Negative <- R_Negative + r[m1 + 1 + t] # r denotes the counting number from -m1 to m2,
    K_Negative <- K_Negative + k[t] # k denotes the counting number from 1 to m2,

    ### Update Parameters for updating Spearman footrule of different imputations in X givn imputations in Y
    # Update Parameters in U
    ranks_Y_U <- ranks_Y_U + 1

    # Update Parameters in Phi
    if( (b > 0) & (a <= 0) ){
      S_Positive <- S_Positive + 1
    }

    if( (-b >= -m1) & (-b <= -1) ){
      # s[as.character(-b)] <- s[as.character(-b)] - 1
      s[m1 + 1 + (-b)] <- s[m1 + 1 + (-b)] - 1    # s denotes the counting number from -m1 to -1,
    }

    if( (-a >= -m1) & (-a <= -1)){
      # s[as.character(-a)] <- s[as.character(-a)] + 1
      s[m1 + 1 + (-a)] <- s[m1 + 1 + (-a)] + 1    # s denotes the counting number from -m1 to -1,
    }

    ## Update Parameters in [n] \ (U \cup Phi)
    # R_Positive <- R_Positive - r[as.character(t-1)]
      R_Positive <- R_Positive - r[m1 + 1 + (t-1)] # r denotes the counting number from -m1 to m2,

    # r_star is the counting number such that r_star[as.character(i)] = sum_{i \in [n] (U \cup Phi)} d_i = i
    r_star <- r
    for (i in seq(-m1,-1)) {
     # r_star[as.character(i)] <- r_star[as.character(i-1)]    ############ Warning: this is incorrect, check here more carefully
     # r_star[as.character(i)] <- r[as.character(i+t)]
      r_star[m1 + 1 + i] <- r[m1 + 1 + i + t]  # r_star denotes the counting number from -m1 to -1,
    }

    ## Calculate the maximum Spearman footrule given imputations in Y
    r_s <- rep(0, m1)
    # names(r_s) <- as.character(seq(-m1,-1))
    #
    # for (i in seq(-m1, -1)) {
    #   r_s[as.character(i)] <- r_star[as.character(i)] + s[as.character(i)]
    # }

    for (i in seq(1, m1)) {
      r_s[i] <- r_star[i] + s[i] # r_star denotes the counting number from -m1 to -1,  # s denotes the counting number from -m1 to -1,
      # r_s denotes the counting number from -m1 to -1
    }

    R_S_Positive <- R_Positive + S_Positive

    Res <- UpperBoundCaseOne(Initial_D = D[t+1,1], ranks_Y_U = ranks_Y_U, R_Positive = R_S_Positive, r = r_s, n = n, m3 = m3, t3 = t3, fullres= TRUE )

    D[t+1,] <- Res

  }

  # Return the maximum possible Spearman footrule distance
  if(fullres == TRUE){
    return(D)
  }else{
    return(max(D))
  }
}




