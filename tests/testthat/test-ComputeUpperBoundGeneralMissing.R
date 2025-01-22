test_that("gives the upper bound as given by
          enumerating algorithm: missing case I", {

            X <- c(1, 5, NA, 3, 2, 4, NA,6)
            Y <- c(2, 7, 6, 1, 3, 5, 4, 8)
            # enumerating
            Res <- ElaborateSpearmanFootrule(X,Y)
            enumerating_max = max(Res$AllPossibleSPearmanFootrule)
            expect_equal(ComputeUpperBoundGeneralMissing(X, Y), enumerating_max)
})


test_that("gives the upper bound as given by
          enumerating algorithm: missing case II", {

            X <- c(1, 5, NA, 3, 2, 4, NA, 6)
            Y <- c(2, NA, 6, 1, 3, 5, 4, NA)
            # enumerating
            Res <- ElaborateSpearmanFootrule(X,Y)
            enumerating_max = max(Res$AllPossibleSPearmanFootrule)
            expect_equal(ComputeUpperBoundGeneralMissing(X, Y), enumerating_max)
})


test_that("gives the upper bound as given by
          enumerating algorithm: missing case III", {
            
            X <- c(1, 5, NA, 3, 2, 4, NA, 6)
            Y <- c(2, 7, NA, 1, 3, 5, NA, 8)
            # enumerating
            Res <- ElaborateSpearmanFootrule(X,Y)
            enumerating_max = max(Res$AllPossibleSPearmanFootrule)
            expect_equal(ComputeUpperBoundGeneralMissing(X, Y), enumerating_max)
})


test_that("gives the upper bound as given by
          enumerating algorithm: general missing case", {
            X <- c(1, 5, NA, 3, 2, 4, NA, 6)
            Y <- c(2, NA, NA, 1, 3, 5, 4, 6)
            # enumerating
            Res <- ElaborateSpearmanFootrule(X,Y)
            enumerating_max = max(Res$AllPossibleSPearmanFootrule)
            expect_equal(ComputeUpperBoundGeneralMissing(X, Y), enumerating_max)
          })