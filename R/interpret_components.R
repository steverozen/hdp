#' Extract components and exposures from multiple posterior sample chains.
#'
#' This function returns components in which we have high confidence, moderate confidence, and little confidence
#' (likely noise).
#'
#' @param multi.chains.retval A list of objects returned from \code{\link[mSigHdp]{CombineChainsAndExtractSigs}}.
#'
#' @param confident.prop Components found in >= \code{confident.prop} of
#'   posterior samples are high confidence components.
#'
#' @param noise.prop Components found in < \code{noise.prop} of posterior samples are considered noise components.
#'
#' @param verbose if TRUE, generate progress messages.
#'
#' @return In the information that follows, a "component" is
#'  the union of multiple raw clusters of mutations (in the case of mutational
#'  signature analysis). Invisibly, a list with the following elements: \describe{
#'
#' \item{high_confident_components}{A data frame containing the
#'    components found in >= \code{confident.prop} of posterior samples. Each column is
#'    a component; in the case of mutational signatures the rows are mutation types.}
#'
#' \item{high_confident_components_post_number}{A data frame
#'   in which the first column contains the index of a column in \code{high_confident_components}
#'  and the second column contains the number of posterior samples that
#'  contributed to that component.}
#'
#' \item{high_confident_components_cdc}{A matrix in which each row
#'  corresponds to one of the Dirichlet processes, and each column
#'  corresponds to one component in \code{high_confident_components}.
#'  In the case of mutational signature analysis, most of the columns
#'  correspond to an input biological sample (e.g. individual tumor).}
#'
#' \item{moderate_components}{Analogous to \code{high_confident_compents} except for
#'  components with constituent raw clusters found in >= \code{noise.prop} and
#'   < \code{confident.prop} posterior samples.}
#'
#' \item{moderate_components_post_number}{Analogous to \code{high_confident_components_post_number}.}
#'
#' \item{moderate_components_cdc}{Analogous to \code{high_confident_components_cdc}.}
#'
#' \item{noise_components}{Analogous to \code{high_confident_compents} except for
#'  components with constituent raw clusters found in  < \code{noise.prop} posterior samples.}
#'
#' \item{noise_components_post_number}{Analogous to \code{high_confident_components_post_number}.}
#'
#' \item{noise_components_cdc}{Analogous to \code{high_confident_components_cdc}.}
#'
#' }
#'
#' @export

interpret_components <- function(multi.chains.retval,
                                 confident.prop = 0.90,
                                 noise.prop = 0.50,
                                 verbose=T) {
  if (verbose) message("extracting components ", Sys.time())

  components_category_counts <- multi.chains.retval$components
  # A data frame in which each column represents a component: ie a cluster of
  # mutations created from the union of some raw components that we think
  # represent the same abstract cluster. In the application to mutational
  # signature extraction, each row is a mutation type; in the application to
  # probabilistic topic modeling, each row is a word. Each cell contains the
  # number of mutations of a given mutation type in a given cluster.

  components_post_number <- multi.chains.retval$components.post.samples
  # A data frame with with 2 columns. The first column contains the index of a
  # column in components_category_counts, the second column contains the number of posterior
  # samples in which the raw components contributing the component appeared.


  nsamp <-  multi.chains.retval$nsamp
  components_cdc <- multi.chains.retval$components.cdc

  components_category_counts <- components_category_counts[,order(components_post_number[,2],decreasing=T),drop=F]
  components_cdc <- components_cdc[,order(components_post_number[,2],decreasing=T),drop=F]
  components_post_number <- components_post_number[order(components_post_number[,2],decreasing=T),]

  #the components with more than confident.prop nsamples are
  #selected as components with high confidence
  high_confident_components <- components_category_counts[,which(components_post_number[,2]>=(confident.prop*nsamp))]
  high_confident_components_post_number <- components_post_number[which(components_post_number[,2]>=(confident.prop*nsamp)),]
  high_confident_components_cdc <- components_cdc[,which(components_post_number[,2]>=(confident.prop*nsamp))]
  #the components with more than noise.prop nsamples but less
  #than confident.prop samples are selected as components with
  #moderate confidence
  moderate_components <- components_category_counts[,intersect(which(components_post_number[,2]>=(noise.prop*nsamp)),which(components_post_number[,2]<(confident.prop*nsamp)))]
  moderate_components_post_number <- components_post_number[intersect(which(components_post_number[,2]>=(noise.prop*nsamp)),which(components_post_number[,2]<(confident.prop*nsamp))),]
  moderate_components_cdc <- components_cdc[,intersect(which(components_post_number[,2]>=(noise.prop*nsamp)),which(components_post_number[,2]<(confident.prop*nsamp)))]

  #the components with less than noise.prop nsamples
  #are selected as noise components
  noise_components <- components_category_counts[,which(components_post_number[,2]<(noise.prop*nsamp))]
  noise_components_post_number <- components_post_number[which(components_post_number[,2]<(noise.prop*nsamp)),]
  noise_components_cdc <- components_cdc[,which(components_post_number[,2]<(noise.prop*nsamp))]

  #noise components also include components that only occur in one posterior
  #sample of each chain
  noise_components <- cbind(noise_components,multi.chains.retval$each.chain.noise.components)

  noise_components_cdc <- cbind(noise_components_cdc,multi.chains.retval$each.chain.noise.cdc)




  return(invisible(list(high_confident_components               = high_confident_components,
                        high_confident_components_post_number   = high_confident_components_post_number,
                        high_confident_components_cdc           = high_confident_components_cdc,
                        moderate_components                     = moderate_components,
                        moderate_components_post_number         = moderate_components_post_number,
                        moderate_components_cdc                 = moderate_components_cdc,
                        noise_components                        = noise_components,
                        noise_components_post_number            = noise_components_post_number,
                        noise_components_cdc                    = noise_components_cdc
                        )))

}

