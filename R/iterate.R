# posterior MCMC sampling for hdp
# returns an R list with two elements:
# the updated hdp state as a list, and
# the vector of likelihoods for these iterations

#' @useDynLib hdpx hdpMultinomial_iterate
iterate <- function(hdplist, numiter, cpiter, verbosity){
  out <- .Call(hdpMultinomial_iterate, hdplist, numiter, cpiter,
               dolik=1, verbosity, PACKAGE="hdp")
  return(out)
}
