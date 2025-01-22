# Compute lower bounds of Spearman's footrule under general missing case
#source('ComputeLowerBoundCaseOne.R')
#source('ComputeLowerBoundCaseTwo.R')
#source('ComputeLowerBoundCaseThree.R')

ComputeLowerBoundGeneralMissing <- function(X, Y) {
  #X, Y: potentially unobserved data

  n <- length(X)

  # Identify missing values in X and Y
  U <- which(is.na(X) & !is.na(Y))
  Phi <- which(is.na(Y) & !is.na(X))
  Psi <- which(is.na(X) & is.na(Y))

  # Number of missing samples
  m1 <- length(U)
  m2 <- length(Phi)
  m3 <- length(Psi)

  ## No missing data
  if((m1 + m2 + m3 == 0)){
    return(list(Minfootrule = sum(abs(rank(X)-rank(Y)))))
  }

  ## Missing case I
  if( (m3 == 0) & (m1 != 0) & (m2 == 0) ){
    Res <- ComputeLowerBoundCaseOne(X,Y)
    return(Res)
  }
  if( (m3 == 0) & (m1 == 0) & (m2 != 0) ){
    Res <- ComputeLowerBoundCaseOne(Y,X)
    return(Res)
  }

  ## Missing case III
  if( (m3 != 0) & (m1 == 0) & (m2 == 0) ){
    Res <- ComputeLowerBoundCaseThree(X,Y)
    return(Res)
  }

  ## Missing case II and general missing cases
    Res <- ComputeLowerBoundCaseTwo(X[setdiff(seq(1,n), Psi)],Y[setdiff(seq(1,n), Psi)])
    return(Res)
}


