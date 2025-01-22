# Function to reorder X and Y vectors
SortingAlg <- function(X, Y) {

  n <- length(X)
  # Identify missing values in X and Y
  U <- which(is.na(X) & !is.na(Y))
  Phi <- which(is.na(Y) & !is.na(X))
  Psi <- which(is.na(X) & is.na(Y))

  # Number of missing samples
  m1 <- length(U)
  m2 <- length(Phi)
  m3 <- length(Psi)

  # Ensure no overlap in missing indexes
  # Sort Y[U] in ascending order
  Y_U <- sort(Y[U], decreasing = FALSE)

  # Sort X[Phi] in ascending order
  X_Phi <- sort(X[Phi], decreasing = FALSE)

  # Sort X and Y based on sorted indices
  SortedX <- c(rep(NA, m1), X_Phi, X[setdiff(1:n, c(U, Phi, Psi))], rep(NA, m3))
  SortedY <- c(Y_U, rep(NA, m2), Y[setdiff(1:n, c(U, Phi, Psi))], rep(NA,m3))

  return(list(SortedX = SortedX, SortedY = SortedY))
}

# Example usage
#X <- c(1, 2, 5, 4, 3, 6) # X vector with missing values as NA
#Y <- c(5, 4, NA, NA, 1, 3) # Y vector with missing values as NA
#result <- SortingAlg(X, Y)

