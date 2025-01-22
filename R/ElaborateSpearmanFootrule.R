### compute all possible values of Spearman's footrule by elaborating
# Note: be aware of the computational time

# Function to compute factorial
factorial <- function(n) {
  if(n == 0) {
    return(1)
  } else {
    return(n * factorial(n-1))
  }
}

# Function to compute Spearman footrule distance
SpearmanFootrule <- function(X, Y) {
  return(sum(abs(rank(X) - rank(Y))))
}

# Function to generate all combination of any n distinct number from vector X
permutations <- function(vector, m) {
  # Use the permutations function from gtools package
  perms <- gtools::permutations(n = length(vector), r = m, v = vector, repeats.allowed = FALSE)
  return(perms)
}


# Function to elaborate all possible values of Spearman footrule
ElaborateSpearmanFootrule <- function(X, Y) {

  n <- length(X)

  missingIndicesX <- which(is.na(X))
  missingIndicesY <- which(is.na(Y))

  observeIndicesX <- setdiff(seq(1,n),missingIndicesX)
  observeIndicesY <- setdiff(seq(1,n),missingIndicesY)
  # Compute ranks
  ranks_X <- rank(X, na.last = 'keep')
  ranks_Y <- rank(Y, na.last = 'keep')

  mx <- length(missingIndicesX)
  my <- length(missingIndicesY)

  # The number of possibilities of imputed X and Y
  Number_Permutation_X <- factorial(n) / factorial(n-mx)
  Number_Permutation_Y <- factorial(n) / factorial(n-my)

  # All possible ranks of imputed missing samples in X and Y
  if(mx != 0){
    AllPossibleImputed_X <- permutations(seq(1,n), mx)
  }else{
    AllPossibleImputed_X <- matrix(data = NA, nrow = 1, ncol = 1)  ## This could be anything
  }

  if(my != 0){
    AllPossibleImputed_Y <- permutations(seq(1,n), my)
  }else{
    AllPossibleImputed_Y <- matrix(data = NA, nrow = 1, ncol = 1)  ## This could be anything
  }

  # Initialize D
  D <- 0
  Count <- 1
  for(i in seq(1,Number_Permutation_X)){

    for (j in seq(1,Number_Permutation_Y)) {

      ## Impute X
      Imputed_X <- ranks_X
      Imputed_X[missingIndicesX] <- AllPossibleImputed_X[i,]
      Imputed_X[observeIndicesX] <- setdiff(seq(1,n), AllPossibleImputed_X[i,])[ranks_X[observeIndicesX]]

      ## Impute Y
      Imputed_Y <- ranks_Y
      Imputed_Y[missingIndicesY] <- AllPossibleImputed_Y[j,]
      Imputed_Y[observeIndicesY] <- setdiff(seq(1,n), AllPossibleImputed_Y[j,])[ranks_Y[observeIndicesY]]

      ## Compute Spearman footrule
      D[Count] <- SpearmanFootrule(Imputed_X, Imputed_Y)
      Count <- Count + 1

    }
  }

  return(list(AllPossibleSPearmanFootrule = D))
}

# # Example usage with small vectors and a small number of missing values
# X <- c(NA, 2, 3, NA, 5)
# Y <- c(1, 4, 2, 6, 3)
# footruleValues <- ElaborateSpearmanFootrule(X, Y)

