pkgname <- "bosfr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "bosfr-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('bosfr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("boundsKendall")
### * boundsKendall

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: boundsKendall
### Title: Bounds of Kendall's tau in the Presence of Missing Data
### Aliases: boundsKendall

### ** Examples

### compute bounds of Kendall's tau between incomplete ranked lists
X <- c(1, 2, NA, 4, 3)
Y <- c(3, NA, 4, 2, 1)
boundsKendall(X, Y)

### compute bounds of Kendall's tau between incomplete vectors of distinct data
X <- c(1.3, 2.6, NA, 4.2, 3.5)
Y <- c(5.5, NA, 6.5, 2.6, 1.1)
boundsKendall(X, Y)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("boundsKendall", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("boundsSFR")
### * boundsSFR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: boundsSFR
### Title: Exact bounds of Spearman's footrule in the Presence of Missing
###   Data
### Aliases: boundsSFR

### ** Examples

### compute exact bounds of Spearman's footrule between incomplete ranked lists
X <- c(1, 2, NA, 4, 3)
Y <- c(3, NA, 4, 2, 1)
boundsSFR(X, Y, pval=FALSE)

### compute exact bounds of Spearman's footrule between incomplete vectors of distinct data,
### and perform independence test
X <- c(1.3, 2.6, NA, 4.2, 3.5)
Y <- c(5.5, NA, 6.5, 2.6, 1.1)
boundsSFR(X, Y, pval=TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("boundsSFR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
