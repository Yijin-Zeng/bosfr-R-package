# Compute upper bounds of Spearman footrule given:
#  Initial_D, R_Positive, r, m3, n

UpperBoundCaseThree <- function(Initial_D, R_Positive, r, m3, n, fullres = FALSE){

  # Initialize D
  D <- rep(0, m3 + 1)

  # Calculate initial Spearman footrule distance
  D[1] <- Initial_D

  # r_MinusOne <- r[as.character(-1)]
  r_MinusOne <- r[2*m3 + 2 - 1]  # r denotes the counting number from (-2*m3 - 1) to -1
  # Iteratively calculate upper bound for each t
  for (t in 1:m3) {
    a <- t - (n + 1 - t)      # the differences of ranks in n - t + 1 after altering
    b <- (n - m3 + t) - (m3 - t + 1) # the differences of ranks in n - t + 1 before altering

    D[t+1] <- D[t] + abs(a) - abs(b) + 2*R_Positive - 2*((n - m3) - R_Positive) + 2*r_MinusOne

    # R_Positive <- R_Positive + as.numeric(r[as.character(-2*t + 1)]) + as.numeric(r[as.character(-2*t)])
    R_Positive <- R_Positive + r[2*m3 + 2 -2*t + 1] + r[2*m3 + 2 -2*t] # r denotes the counting number from (-2*m3 - 1) to -1

    #r_MinusOne <- r[as.character( (-2*t -1) )]
    r_MinusOne <- r[2*m3 + 2 -2*t -1]
  }

  if(fullres == TRUE){
    return(D)
  }else{
    return(max(D))
  }

}
