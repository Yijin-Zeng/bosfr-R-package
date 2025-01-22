#' @md
#' @export
#'
#' @title Bounds of Kendall's tau in the Presence of Missing Data
#'
#' @description Computes bounds of Kendall's tau in the
#' presence of missing data. Suitable only for univariate distinct data
#' where no ties is allowed.
#'
#' @usage boundsKendall(X, Y)
#'
#' @param X,Y Numeric vectors of data values with potential missing data.
#' No ties in the data is allowed. Inf and -Inf values will be omitted.
#'
#' @details \code{boundsKendall()} computes bounds of Kendall's tau
#' for partially observed univariate, distinct data. The bounds are computed
#' by first calculating the bounds of Spearman's footrule (\cite{Zeng et al., 2025}), and then applying
#' the combinatorial inequality between Kendall's tau and Spearman's footrule
#' (\cite{Kendall, 1948}). See \cite{Zeng et al., 2025} for more details.
#'
#' Let \eqn{X = (x_1, \ldots, x_n)} and \eqn{Y = (y_1, \ldots, y_n)} be
#' two vectors of univariate, distinct data.
#' Kendall's tau is defined as the number of discordant pairs between \eqn{X} and \eqn{Y}:
#' \deqn{\tau(X,Y) = \sum\limits_{i < j} \{I(x_i < x_j)I(y_i > y_j) + I(x_i > x_j)I(y_i < y_j)\}.}
#' Scaled Kendall's tau \eqn{\tau_{Scale}(X,Y) \in [0,1]} is defined as (\cite{Kendall, 1948}):
#' \deqn{\tau_{Scale}(X,Y) = 1 - 4\tau(X,Y)/(n(n-1)).}
#'
#' @return
#'
#'  \item{bounds}{bounds of Kendall's tau.}
#'
#'  \item{bounds.scaled}{bounds of scaled Kendall's tau.}
#'
#' @references
#' \itemize{
#'  \item Zeng Y., Adams N.M., Bodenham D.A. Exact Bounds of Spearman's footrule in the Presence of Missing Data with Applications to Independence Testing. arXiv preprint arXiv:2501.11696. 2025 Jan 20.
#'  \item Kendall, M.G. (1948) Rank Correlation Methods. Charles Griffin, London.
#'  \item Diaconis, P. and Graham, R.L., 1977. Spearman's footrule as a measure of disarray. Journal of the Royal Statistical Society Series B: Statistical Methodology, 39(2), pp.262-268.
#' }
#'
#' @examples
#' ### compute bounds of Kendall's tau between incomplete ranked lists
#' X <- c(1, 2, NA, 4, 3)
#' Y <- c(3, NA, 4, 2, 1)
#' boundsKendall(X, Y)
#'
#' ### compute bounds of Kendall's tau between incomplete vectors of distinct data
#' X <- c(1.3, 2.6, NA, 4.2, 3.5)
#' Y <- c(5.5, NA, 6.5, 2.6, 1.1)
#' boundsKendall(X, Y)

boundsKendall <- function(X, Y){

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
  upper_bound_footrule <- ComputeUpperBoundGeneralMissing(X,Y)
  lower_bound_footrule <- ComputeLowerBoundGeneralMissing(X,Y)$Minfootrule

  upper_bound = upper_bound_footrule
  lower_bound = lower_bound_footrule/2

  # scaled upper and lower bounds
  n <- length(X)
  scaled_upper_bound <- min(1,1 - 4*lower_bound/(n*(n-1)))
  scaled_lower_bound <- max(-1,1 - 4*upper_bound/(n*(n-1)))

  RES <- list(bounds = c(lower_bound, upper_bound),
              bounds.scaled = c(scaled_lower_bound, scaled_upper_bound))

  return(RES)
}
