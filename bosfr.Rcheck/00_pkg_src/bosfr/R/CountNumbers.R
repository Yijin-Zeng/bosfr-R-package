# Computes the values of samples in X included in N
CountNumbers <- function(X,N){
# X is the vector which we want to count the numbers in it
# N is the vector including the numbers we want count, N should be a, a+1, ..., b

  # Initialize
  number <- length(N)
  Res <- rep(0, number)
  # names(Res) <- as.character(N)

  a <- N[1]
  b <- N[number]
  X_filtered <- X[X >= a & X <= b]  ## Only consider samples between a and b

  for(i in X_filtered){
    # Res[as.character(i)] <- Res[as.character(i)] + 1
    Res[i - a + 1] <- Res[i - a + 1] + 1

  }

  return(list(Counted_Numbers = Res))
}

# ## Examples
# X <- c(0,1,2,3,-5,-2,-2,1)
# N <- c(-3,-2)
# CountNumber(X,N)


