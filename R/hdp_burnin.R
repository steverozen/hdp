#' Burnin of posterior sampling chain across activated DPs.
#'
#' Run a Gibbs sampler over the activated nodes of a hierarchical Dirichlet process.
#' Each iteration re-assigns the cluster allocation of every data item.
#'
#' @param hdp An \code{\link[hdpx]{hdpState-class}} object.
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
#' \item{hdplist}{A list represetnation of an \code{\link[hdpx]{hdpState-class}}.}
#'
#' \item{liklihood}{A numeric vector with the likelihood at each iteration.}
#'
#' }
#'
#' @export
#'
#' @examples
#' my_hdp <- hdp_init(ppindex=0, cpindex=1, hh=rep(1, 6), alphaa=rep(1, 3), alphab=rep(2, 3))
#' my_hdp <- hdp_adddp(my_hdp, 2, 1, 2)
#' my_hdp <- hdp_adddp(my_hdp, 10, c(rep(2, 5), rep(3, 5)), 3)
#' my_hdp <- hdp_setdata(my_hdp, 4:13, example_data_hdp)
#' my_hdp <- dp_activate(my_hdp, 1:13, 2)
#' my_hdp_chain <- hdp_burnin(my_hdp, burnin = 100, cpiter = 100)
hdp_burnin <- function(hdp,
                       burnin,
                       cpiter=1,
                       verbosity=0){

  # input checks
  if (class(hdp) != "hdpState") stop("hdp must have class hdpState")
  if (!validObject(hdp)) stop("input hdp is not a valid hdpState object")
  for (arg in c("burnin", "cpiter")) {
    x <- get(arg)
    if (x < 1 | x %% 1 != 0) stop(paste(arg, "must be a positive integer"))
  }
  if (verbosity < 0 |
        verbosity > 4 |
        verbosity %% 1 != 0) stop("verbosity must be integer from 0--4")

  # translate hdp hdpState (S4 class) to plain list so C code can parse
  hdplist <- as.list(hdp)

  burnin.time <- system.time(
    # Call the C code to do the Gibbs sampling iterations
    output <- iterate(hdplist, burnin, cpiter, verbosity)
  )
  message("hdp_burnin time: ")
  for (xn in names(burnin.time)) {
    message(" ", xn, " ", burnin.time[[xn]])
  }

  return(invisible(list(hdplist     = output[[1]],
                       likelihood = output[[2]])))

}
