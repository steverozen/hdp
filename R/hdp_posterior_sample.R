#' Posterior sampling chain across activated DPs.
#'
#' Run a Gibbs sampler over the burnin chains from \code{\link{hdp_burnin}}.
#' Each iteration re-assigns the cluster allocation of every data item.
#' Run \code{burnin} iterations, and then collect \code{n} samples from the chain
#' with \code{space} iterations between each collected sample. To plot output,
#' see \code{\link{plot_lik}}, \code{\link{plot_numcluster}}, and
#' \code{\link{plot_data_assigned}}. Can collect multiple
#' independent HDP sampling chains in a hdpSampleMulti object via \code{\link{hdp_multi_chain}}.
#' Components are extracted via \code{\link{hdp_extract_components}}.
#'
#' @param burnin.output An S4 object from \code{\link{hdp_burnin}}.
#' @param n The number of posterior samples to collect.
#' @param space The number of iterations between collected samples.
#' @param cpiter The number of iterations of concentration parameter sampling to perform after each iteration.
#' @param seed The (integer) seed that can be set to reproduce output. Default is a
#'  random seed from 1 -- 10^7, reported in the output.
#' @param verbosity Verbosity of debugging statements.
#'  0 (least verbose) -- 4 (most verbose). 0 highly recommended - only change for debugging small examples.
#' @return A hdpSampleChain object with the salient information from each
#'  posterior sample. See \code{\link{hdpSampleChain-class}}
#' @seealso \code{\link{hdp_multi_chain}}, \code{\link{hdp_extract_components}},
#'  \code{\link{cull_posterior_samples}}, \code{\link{plot_lik}}, \code{\link{plot_numcluster}},
#'  \code{\link{plot_data_assigned}}
#' @importClassesFrom Matrix dgCMatrix
#' @export

hdp_posterior_sample <- function(burnin.output,
                                 n,
                                 space,
                                 cpiter=1,
                                 seed=sample(1:10^7, 1),
                                 verbosity=0){

  set.seed(seed) ##set.seed in the function because of parallelization

  # input checks
  ## check burnin.output

  hdplist <- burnin.output$hdplist

  if (cpiter < 1 | cpiter %% 1 != 0) stop("cpiter must be a positive integer")

  if (verbosity < 0 |
      verbosity > 4 |
      verbosity %% 1 != 0) stop("verbosity must be integer from 0--4")

  # initialise list for posterior sample output
  starttime <- Sys.time()
  curriter <- 0 #keep track of time
  sample  <- list()

  if(!exists("all.lik",where=burnin.output)){
    burnin.output$all.lik <- burnin.output$likelihood
  } ##maybe all.lik is not available for every burnin? not sure. but safe to check

  all.lik <- burnin.output$all.lik
  burnin <- length(burnin.output$all.lik) #pass burnin to the final state
  totiter <- n*space
  # collect n posterior samples
  for (samp in 1:n){

    output <- iterate(hdplist, space, cpiter, verbosity)
    hdplist <- output[[1]]
    all.lik <- c(all.lik,output[[2]])

    sample <- c(sample,
                list(hdp_getstate(hdplist)))


    #report time every 10 samples if > 1 min has passed
    tracktime <- Sys.time()
    curriter <- curriter + space

    if (samp %% 10 == 0){
      elapsedtime <- (tracktime - starttime)/60
      print(sprintf("time %1.1f ETC %1.1f mins",
                    elapsedtime, elapsedtime / curriter * totiter))
    }
  }


  numclass <- sapply(sample, function(x) x$numclass)
  classqq <- lapply(sample, function(x) x$classqq)
  classnd <- lapply(sample, function(x) as(x$classnd, "dgCMatrix"))
  alpha <- t(sapply(sample, function(x) x$alpha))

  # if only one conparam, then alpha can have wrong dims (vector not matrix)
  # Need to check. I don't think there is any condition has only one conparam. added by Mo
  if (dim(alpha)[1]==1 & n > 1) {
    alpha <- matrix(alpha, ncol=1)
  }

  #translate hdplist back to HDPObject class
  hdp <- as.hdpState(hdplist)
  remove(hdplist)

  ans <- new("hdpSampleChain",
             seed = as.integer(seed),
             settings = list(burnin=burnin,
                             n=n,
                             space=space,
                             cpiter=cpiter),
             hdp = hdp,
             lik = all.lik,
             numcluster = numclass,
             cp_values = alpha,
             clust_categ_counts = classqq,
             clust_dp_counts = classnd,
             numcomp = as.integer(NULL),
             prop.ex = as.numeric(NULL),
             comp_cos_merge = as.numeric(NULL),
             comp_categ_counts = list(),
             comp_dp_counts = list(),
             comp_categ_distn = list(),
             comp_dp_distn = list())

  # check validity and return
  if (!validObject(ans)) warning("Not a valid hdpSampleChain object.")
  return(ans)
}
