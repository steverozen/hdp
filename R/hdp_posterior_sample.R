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
#' @param post.input An S4 object from \code{\link{hdp_burnin}}.
#' @param post.n The number of posterior samples to collect.
#' @param post.space The number of iterations between collected samples.
#' @param post.cpiter The number of iterations of concentration parameter sampling to perform after each iteration.
#' @param seed The (integer) seed that can be set to reproduce output. Default is a
#'  random seed from 1 -- 10^7, reported in the output.
#' @param post.verbosity Verbosity of debugging statements.
#'  0 (least verbose) -- 4 (most verbose). 0 highly recommended - only change for debugging small examples.
#' @param checkpoint If \code{TRUE}, a checkpoint will be saved for every 10 posterior samples
#' @return A hdpSampleChain object with the salient information from each
#'  posterior sample. See \code{\link{hdpSampleChain-class}}
#' @seealso \code{\link{hdp_multi_chain}}, \code{\link{hdp_extract_components}},
#'  \code{\link{cull_posterior_samples}}, \code{\link{plot_lik}}, \code{\link{plot_numcluster}},
#'  \code{\link{plot_data_assigned}}
#' @importClassesFrom Matrix dgCMatrix
#' @export

hdp_posterior_sample <- function(post.input,
                                 post.n,
                                 post.space,
                                 post.cpiter=1,
                                 seed=sample(1:10^7, 1),
                                 post.verbosity=0,
                                 checkpoint=F){

  set.seed(seed) ##set.seed in the function when running parallel

  ## check the class of input from hdp_burnin or hdp_posterior_sample

  if(class(post.input)[1] == "hdpSampleChain"){ #extend Gibbs sampling

    message("Extend Gibbs Sampling")

    post.input <- as.list(post.input)

    #pick up from the end of last gibbs sampling
    sampling.input <- list(hdplist = as.list(post.input$hdp),
                          likelihood = post.input$lik)

  }else if(class(post.input) == "list"){
      sampling.input <- post.input
      message("Gibbs Sampling after Burn-in Iteration")

    }else{
    message("Input is not a hdpSampleChain or a list from hdp_burnin")
    stop()
    }

  remove(post.input)

  # input checks
  ## check sampling.input

  hdplist <- sampling.input$hdplist

  if (post.cpiter < 1 | post.cpiter %% 1 != 0) stop("post.cpiter must be a positive integer")

  if (post.verbosity < 0 |
      post.verbosity > 4 |
      post.verbosity %% 1 != 0) stop("post.verbosity must be integer from 0--4")

  # initialise list for posterior sample output
  starttime <- Sys.time()
  curriter <- 0 #keep track of time
  sample  <- list()

  if(!exists("all.lik",where = sampling.input)){
    sampling.input$all.lik <- sampling.input$likelihood
  } ##all.lik is not available every time

  all.lik <- sampling.input$all.lik
  burnin <- length(sampling.input$all.lik) #pass burnin to the final state
  totiter <- post.n * post.space

  # function to return difference in time in minute units
  mindifftime <- function(t1, t2){
    as.numeric(t2-t1, units="mins")
  }

  prevtime <- Sys.time()
  # collect post.n posterior samples
  for (samp in 1:post.n){

    output <- iterate(hdplist, post.space, post.cpiter, post.verbosity)
    hdplist <- output[[1]]
    all.lik <- c(all.lik,output[[2]])

    sample <- c(sample,
                list(hdp_getstate(hdplist)))


    #report time every 10 samples if > 1 min has passed
    tracktime <- Sys.time()
    curriter <- curriter + post.space
    if (mindifftime(prevtime, tracktime) > 1 & samp %% 10 == 0){
      elapsedtime <- mindifftime(starttime, tracktime)
      print(sprintf("Current sampling: %1.1f mins;Estimated Time of Completion: %1.1f mins",
                    elapsedtime, elapsedtime / curriter * totiter))
      prevtime <- tracktime
    }
    if(checkpoint){
      if(samp %% 10 == 0){
        message("Checkpoint after collecting ",samp," posterior samples")
        numclass <- sapply(sample, function(x) x$numclass)
        classqq <- lapply(sample, function(x) x$classqq)
        classnd <- lapply(sample, function(x) as(x$classnd, "dgCMatrix"))
        alpha <- t(sapply(sample, function(x) x$alpha))

        hdp <- as.hdpState(hdplist)

        ans <- new("hdpSampleChain",
                   seed = as.integer(seed),

                   ##changing names here cause a check error in hdpx::new
                   settings = list(burnin       = burnin,
                                   n            = samp,
                                   space        = post.space,
                                   cpiter       = post.cpiter),
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
        save(ans,
             file = paste0("posterior.checkpoint.", seed, ".Rdata"))

      }
    }
  }

  ## these four can be wrapped up?
  numclass <- sapply(sample, function(x) x$numclass)
  classqq <- lapply(sample, function(x) x$classqq)
  classnd <- lapply(sample, function(x) as(x$classnd, "dgCMatrix"))
  alpha <- t(sapply(sample, function(x) x$alpha))

  # if only one conparam, then alpha can have wrong dims (vector not matrix)
  # Need to check. I don't think there is any condition has only one conparam. added by Mo
  if (dim(alpha)[1]==1 & post.n > 1) {
    alpha <- matrix(alpha, ncol=1)
  }

  #translate hdplist back to HDPObject class
  hdp <- as.hdpState(hdplist)
  remove(hdplist)

  ans <- new("hdpSampleChain",
             seed = as.integer(seed),

             ##changing names here cause a check error in hdpx::new
             settings = list(burnin       = burnin,
                             n            = post.n,
                             space        = post.space,
                             cpiter       = post.cpiter),
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
