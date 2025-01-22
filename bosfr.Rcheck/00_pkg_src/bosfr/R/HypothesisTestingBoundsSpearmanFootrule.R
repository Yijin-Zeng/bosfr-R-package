### compute the minimum and maximum p-values given the bounds of Speaman's footrule

HypothesisTestingBoundsSpearmanFootrule <- function(min_footrule, max_footrule,n) {
  mu <- 1/3 * n * n
  sigma_square <- 2/45 * n * n * n
  p_1 <- stats::pnorm(min_footrule, mean = mu, sd = sigma_square^(1/2))
  p_2 <- stats::pnorm(max_footrule, mean = mu, sd = sigma_square^(1/2))
  p_1 <- 2*min(p_1, 1-p_1)
  p_2 <- 2*min(p_2, 1-p_2)
  p_min <- min(p_1, p_2)
  p_max <- max(p_1, p_2)
  if( (min_footrule - mu)*(max_footrule - mu) < 0){
    p_max <- 1
  }
  return(c(p_min, p_max))

}
