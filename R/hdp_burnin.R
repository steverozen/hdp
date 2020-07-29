#' Burnin of posterior sampling chain across activated DPs.
#'
#' Run a Gibbs sampler over the activated nodes of a hierarchical Dirichlet process.
#' Each iteration re-assigns the cluster allocation of every data item.
#'
#' @param hdp An \code{\link[hdpx]{hdpState-class}} object or a list representation
#'   of this. The list representation contains elements that correspond to
#'   slots in an\code{\link[hdpx]{hdpState-class}} object.
#'
#' @param burnin The number of burn-in iterations.
#'
#' @param cpiter The number of iterations of concentration parameter sampling
#'  to perform after each iteration.
#'
#' @param verbosity Verbosity of debugging statements.
#'  0 (least verbose) -- 4 (most verbose). 0 highly recommended
#'   - only change for debugging small examples.
#'
#' @return A list with the elements: \describe{
#'
#' \item{hdplist}{A list representation of
#'    an \code{\link[hdpx]{hdpState-class}} object.}
#'
#' \item{likelihood}{A numeric vector with the likelihood at each iteration.}
#'
#' }
#'
#' @export
#'


hdp_burnin <- function(hdp,
                       burnin,
                       cpiter=1,
                       verbosity=0){

  # input checks
  if (class(hdp) == "hdpState") {
    if (!validObject(hdp)) stop("input hdp is not a valid hdpState object")
    # translate hdp hdpState (S4 class) to plain list so C code can parse
    hdplist <- as.list(hdp)
  }  else {
    hdplist <- hdp
  }
  for (arg in c("burnin", "cpiter")) {
    x <- get(arg)
    if (x < 1 | x %% 1 != 0) stop(paste(arg, "must be a positive integer"))
  }
  if (verbosity < 0 |
        verbosity > 4 |
        verbosity %% 1 != 0) stop("verbosity must be integer from 0--4")

  burnin.time <- system.time(
    # Call the C code to do the Gibbs sampling iterations
    output <- iterate(hdplist, burnin, cpiter, verbosity)
  )
  message("hdp_burnin time: ")
  for (xn in names(burnin.time)) {
    message(" ", xn, " ", burnin.time[[xn]])
  }

  return(invisible(list(hdplist   = output[[1]],
                       likelihood = output[[2]])))

}
