#' @md
#' @export
#'
#' @title Exact bounds of Spearman's footrule in the Presence of Missing Data
#'
#' @description Computes exact bounds of Spearman's footrule in the
#' presence of missing data, and performs independence test based on
#' the bounds with controlled Type I error regardless of the values of
#' missing data. Suitable only for univariate distinct data
#' where no ties is allowed.
#'
#' @usage boundsSFR(X, Y, pval = TRUE)
#'
#' @param X Numeric vector of data values with potential missing data.
#' No ties in the data is allowed. Inf and -Inf values will be omitted.
#'
#' @param Y Numeric vector of data values with potential missing data.
#' No ties in the data is allowed. Inf and -Inf values will be omitted.
#'
#' @param pval Boolean for whether to compute the bounds of p-value or not.
#'
#' @details \code{boundsSFR()} computes exact bounds of Spearman's footrule
#' for partially observed univariate, distinct data using the results and
#' algorithms following \cite{Zeng et al., 2025}.
#'
#' Let \eqn{X = (x_1, \ldots, x_n)} and \eqn{Y = (y_1, \ldots, y_n)} be
#' two vectors of univariate, distinct data, and denote the rank of \eqn{x_i}
#' in \eqn{X} as \eqn{R(x_i, X)}, the rank of \eqn{y_i} in \eqn{Y} as
#' \eqn{R(y_i, Y)}.
#' Spearman's footrule is defined as the absolute distance between the ranked
#' values of \eqn{X} and \eqn{Y}:
#' \deqn{D(X,Y) = \sum_{i=1}^{n} |R(x_i, X) - R(y_i, Y)|.}
#' Scaled Spearman's footrule is defined as:
#' \deqn{D_{Scale}(X,Y) = 1 - 3D(X,Y)/(n^2-1).}
#' When \eqn{n} is odd, \eqn{D_{Scale}(X,Y) \in [-0.5,1]}, but when \eqn{n} is
#' even, \eqn{D_{Scale}(X,Y) \in [-0.5\{1+3/(n^2-1)\},1]} (\cite{Kendall, 1948}).
#'
#' The p-value of the independence test using Spearman's footrule, denoted
#' as \eqn{p}, is computed using the normality approximation result in \cite{Diaconis, P., & Graham, R. L. (1977)}.
#' If \code{pval = TRUE}, bounds of the p-value, \eqn{p_{l}, p_{u}} will be
#' computed in the presence of missing data, such that \eqn{p \in [p_{l}, p_{u}]}.
#' The independence test method proposed in \cite{Zeng et al., 2025} returns \eqn{p_{u}} as its p-value.
#' This method controls the Type I error regardless of the values of missing data.
#' See \cite{Zeng et al., 2025} for details.
#'
#' @return
#'
#'  \item{bounds}{exact bounds of Spearman's footrule.}
#'
#'  \item{bounds.scaled}{exact bounds of scaled Spearman's footrule.}
#'
#'  \item{pvalue}{the p-value for the test. (Only present if argument \code{pval = TRUE}.)}
#'
#'  \item{bounds.pvalue}{bounds of the p-value of independence test using Spearman's footrule. (Only present if argument \code{pval = TRUE}.)}
#'
#'
#' @references
#' \itemize{
#'  \item Zeng Y., Adams N.M., Bodenham D.A. Exact Bounds of Spearman's footrule in the Presence of Missing Data with Applications to Independence Testing. arXiv preprint arXiv:2501.11696. 2025 Jan 20.
#'  \item Kendall, M.G. (1948) Rank Correlation Methods. Charles Griffin, London.
#'  \item Diaconis, P. and Graham, R.L., 1977. Spearman's footrule as a measure of disarray. Journal of the Royal Statistical Society Series B: Statistical Methodology, 39(2), pp.262-268.
#' }
#'
#' @examples
#' ### compute exact bounds of Spearman's footrule between incomplete ranked lists
#' X <- c(1, 2, NA, 4, 3)
#' Y <- c(3, NA, 4, 2, 1)
#' boundsSFR(X, Y, pval=FALSE)
#'
#' ### compute exact bounds of Spearman's footrule between incomplete vectors of distinct data,
#' ### and perform independence test
#' X <- c(1.3, 2.6, NA, 4.2, 3.5)
#' Y <- c(5.5, NA, 6.5, 2.6, 1.1)
#' boundsSFR(X, Y, pval=TRUE)

boundsSFR <- function(X, Y, pval=TRUE){

  # Remove all infinite
  X <- X[is.finite(X) | is.na(X)]
  Y <- Y[is.finite(Y) | is.na(Y)]

  # Observed data in X and Y, respectively
  X_prime <- X[!is.na(X)]
  Y_prime <- Y[!is.na(Y)]

  if ((anyDuplicated(X_prime) != 0) | (anyDuplicated(Y_prime) != 0)){
    stop('no Ties is allowed in X or Y')
  }

  # compute upper and lower bounds
  upper_bound <- ComputeUpperBoundGeneralMissing(X,Y)
  lower_bound <- ComputeLowerBoundGeneralMissing(X,Y)$Minfootrule

  # scaled upper and lower bounds
  n <- length(X)
  scaled_upper_bound <- 1-3*lower_bound/(n^2-1)
  scaled_lower_bound <- 1-3*upper_bound/(n^2-1)

  # compute bounds of p-values
  if (pval == TRUE){
    pvalues <- HypothesisTestingBoundsSpearmanFootrule(lower_bound,upper_bound,n)
    min_pvalue <- pvalues[1]
    max_pvalue <- pvalues[2]
    RES <- list(bounds = c(lower_bound, upper_bound),
                bounds.scaled = c(scaled_lower_bound, scaled_upper_bound),
                pvalue = max_pvalue,
                bounds.pvalue = c(min_pvalue, max_pvalue))
  }else{
    RES <- list(bounds = c(lower_bound, upper_bound),
                bounds.scaled = c(scaled_lower_bound, scaled_upper_bound))
  }

  return(RES)
}
