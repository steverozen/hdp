# Caluculate the unsigned Stirling numbers of the first kind up to nn
#
# The return value is the vector of Stirling numbers, each divided by
# the maxium Stirling number in the series. In the comments below
# we will just refer to these as vectors of Stirling numbers.
#
stirling2 <- function(nn){

  # More R notes: parent.frame()returns an environment (termed
  # a "frame" in S). Not sure where parent.frame(2) if this
  # is deep in a call stack. The argument of parent.frame
  # maxes out to .GlobalEnv. This might affect memory usage
  # if different calls are "memo-izing" this function separately.
  # Note that there are no automated tests for package hdp.
  # Perhaps this could be handled better as a closure?

  if (!exists("maxnn", where = parent.frame(2))) {
    assign("maxnn", 1, envir = parent.frame(2))
    assign("allss", list(1), envir = parent.frame(2))
  }

  maxnn.local <- get("maxnn", envir = parent.frame(2))

  #only calculate this if needed
  if (nn > maxnn.local) {

    # allss is a list of all Stirling number series up to maxnn
    allss.local <- get("allss", envir = parent.frame(2))

    allss.local[(length(allss.local) + 1):nn] <- 0 # Initialize each new vector of Stirling number series to a single 0

    for (mm in (maxnn.local + 1):nn) { # For each new Stirling number series

      allss.local[[mm]] <-
        c(allss.local[[mm - 1]] * (mm - 1), 0                    ) +
        c(0,                                allss.local[[mm - 1]])

      mss <- max(allss.local[[mm]])
      allss.local[[mm]] <- allss.local[[mm]] / mss
    }

    assign("maxnn", nn, envir=parent.frame(2))
    assign("allss", allss.local, envir=parent.frame(2))
  }

  assign("nn", nn, envir=parent.frame(2))
  ss <- eval(quote(allss[[nn]]), envir=parent.frame(2))
  # gc()
  return(ss)
}
