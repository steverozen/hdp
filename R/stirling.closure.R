# Caluculate the unsigned Stirling numbers of the first kind up to nn
#
# The return value is the vector of Stirling numbers, each divided by
# the maxium Stirling number in the series. In the comments below
# we will just refer to these as vectors of Stirling numbers.
#
make.stirling <- function(){

  allss <- list(1)

  stir <- function(nn) {
    cat("c", nn, "\n", sep = "")

    len.all <- length(allss)

    for (mm in (len.all + 1):nn) { # For each new Stirling number series, if any

      ss <- allss[[mm - 1]]

      newss <- c(ss * (mm - 1), 0) + c(0, ss)
      newss <- newss / max(newss)

      last.non.0 <- max(which(newss > 0))
      newss <- newss[1:last.non.0]

      allss[[mm]] <<- newss
    }

    retval <- allss[[nn]]
    retval <- c(retval, rep(0, nn - length(retval)))
    stopifnot(length(retval) == nn)
    return(retval)
  }

  return(stir)

}

tests <- function(deb = FALSE) {

  bb <- make.stirling()
  if (deb) debug(bb)
  testthat::expect_equal(bb(3), stirling(3))
  testthat::expect_equal(bb(222), stirling(222))

}
