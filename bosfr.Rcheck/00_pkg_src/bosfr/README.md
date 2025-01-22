
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bosfr (Bounds of Spearman’s Footrule)

<!-- badges: start -->
<!-- badges: end -->

bosfr computes exact bounds of Spearman’s footrule in the presence of
missing data, and performs independence test based on the bounds with
controlled Type I error regardless of the values of missing data. Note
this package is suitable only for univariate distinct data where no ties
is allowed. Bounds of Kendall’s tau is also available with missing data,
provided the bounds of Spearman’s footrule. See Zeng et al., 2025 for
more details.

## Installation

You can install the development version of bosfr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Yijin-Zeng/bosfr-R-package")
```

## Example

There are some basic examples which shows you how to use this package:

``` r
library(bosfr)
### compute exact bounds of Spearman's footrule between incomplete ranked lists
X <- c(1, 2, NA, 4, 3)
Y <- c(3, NA, 4, 2, 1)
boundsSFR(X, Y, pval=FALSE)
#> $bounds
#> [1]  6 12
#> 
#> $bounds.scaled
#> [1] -0.50  0.25
```

``` r

### compute exact bounds of Spearman's footrule between incomplete vectors of distinct data

X <- c(1.3, 2.6, NA, 4.2, 3.5)
Y <- c(5.5, NA, 6.5, 2.6, 1.1)
boundsSFR(X, Y, pval=TRUE)
#> $bounds
#> [1]  6 12
#> 
#> $bounds.scaled
#> [1] -0.50  0.25
#> 
#> $pvalue
#> [1] 1
#> 
#> $bounds.pvalue
#> [1] 0.1197949 1.0000000
```

``` r


### Compute bounds of Kendall's tau
X <- c(1, 2, NA, 4, 3)
Y <- c(3, NA, 4, 2, 1)
boundsKendall(X, Y)
#> $bounds
#> [1]  3 12
#> 
#> $bounds.scaled
#> [1] -1.0  0.4
```

## References

Zeng Y., Adams N.M., Bodenham D.A. Exact Bounds of Spearman’s footrule
in the Presence of Missing Data with Applications to Independence
Testing. arXiv preprint arXiv:2501.11696. 2025 Jan 20.
