#' This function extract signatures from raw clusters discovered by hdp.
#' This returns clusters with high confidence -extracted by more than 90% posterior samples;
#'                            moderate confidence - extracted less than 90% but more than number of post.samples on one chain
#'                            noise - less than number of post.samples on one chain
#' and they are found in how many post.samples in 'stats'
#'
#' ccc_0 and cdc_0 are for diagnostic checking and exposure summary
#'
#' @param x hdpSampleChain or hdpSampleMulti object
#' @param cos.merge Merge components with cosine similarity above this threshold (default 0.90)
#'
#' @return A hdpSampleChain or hdpSampleMulti object updated with component information
#' @aliases extract_sigs_from_clusters
#' @seealso \code{\link{hdp_posterior}}, \code{\link{hdp_multi_chain}},
#'  \code{\link{plot_comp_size}}, \code{\link{plot_comp_distn}},
#'  \code{\link{plot_dp_comp_exposure}}
#' @import clue
#' @export
# @examples
# hdp_extract_components(mut_example_multi)

extract_sigs_from_clusters <-  function(x,
                                        cos.merge = 0.90){
  if (class(x)=="hdpSampleChain") {
    message('Extracting components on single chain.A hdpSampleMulti object is recommended, see ?hdp_multi_chain')
    is_multi <- FALSE
  } else if (class(x)=="hdpSampleMulti") {
    is_multi <- TRUE
  } else {
    stop("x must have class hdpSampleChain or hdpSampleMulti")
  }

  # if (!validObject(x)) stop("x not a valid object")

  if (class(cos.merge) != "numeric" | cos.merge >=1 | cos.merge <= 0) {
    stop("cos.merge must be between 0 and 1")
  }


  if (is_multi) {
    # list of hdpSampleChain objects
    chlist <- x@chains
    nch <- length(chlist)

    # set seed, get final state and number of posterior samples
    set.seed(sampling_seed(chlist[[1]]), kind="Mersenne-Twister", normal.kind="Inversion")
    finalstate <- final_hdpState(chlist[[1]])
    nsamp <- sum(sapply(chlist, function(x) hdp_settings(x)$n))

  }

  # number of categories, DPs,data items at each DP, and frozen priors
  ncat <- numcateg(finalstate) ##number of channel
  ndp <- numdp(finalstate) ##number of dp
  numdata <- sapply(dp(finalstate), numdata) #number of mutations in each sample
  pseudo <- pseudoDP(finalstate)
  rm(finalstate)

  is_prior <- length(pseudo) > 0
  if (is_prior) {
    priorcc <- 1:length(pseudo)
  }


  # Step (1)
  # Make each ccc (clust_categ_counts) and
  # cdc (clust_dp_counts) matrix have the
  # same number of columns
  if(is_multi){
    ccc_0 <- lapply(chlist, function(ch){
      lapply(clust_categ_counts(ch), function(x){
        ans <- cbind(x)
        return(ans[, -ncol(ans)])
      })
    })


    cdc_0 <- lapply(chlist, function(ch){
      lapply(clust_dp_counts(ch), function(x){
        ans <- cbind(x)
        return(ans[, -ncol(ans)])
      })
    })

    # if priors, remove pseudo-counts from ccc_0
    if (is_prior){
      pseudodata <- sapply(dp(final_hdpState(chlist[[1]]))[pseudo],
                           function(x) table(factor(x@datass, levels=1:ncat)))

      ccc_0 <- lapply(ccc_0, function(y) lapply(y, function(x) {
        x[,priorcc] <- x[,priorcc] - pseudodata
        return(x)
      }))
    }
  }


  ########################################################################
  #######merge clusters with high cos.sim in one posterior sample#########
  ########################################################################
  first.merge <- function(ccc,cdc){
    clust_label <- 1:ncol(ccc)
    clust_cos <- lsa::cosine(ccc)
    clust_same <- (clust_cos > 0.95 & lower.tri(clust_cos))
    same <- which(clust_same, arr.ind=TRUE) # merge these columns
    if (length(same)>0){
      for (index in 1:nrow(same)){
        clust_label[same[index, 1]] <- clust_label[same[index, 2]]
      }
    }

    ccc <- merge_cols(ccc,clust_label)
    cdc <- merge_cols(cdc,clust_label)
    return(list(ccc=ccc,cdc=cdc))
  }


  for(i in 1:length(ccc_0)){

    len <- length(ccc_0[[i]])

    test <- mapply(first.merge,ccc_0[[i]],cdc_0[[i]])

    for(j in 1:len){
      ccc_0[[i]][[j]] <- test[[(2*j)-1]]
      cdc_0[[i]][[j]] <- test[[2*j]]
    }



  }
  ########################################################################
  #######merge clusters with high cos.sim in one posterior chain##########
  ########################################################################

  cosmergechain <- function(ccc){
    ccc_unlist <- do.call(cbind,ccc)
    test <- cosCpp(as.matrix(ccc_unlist))
    clust_cos <- test
    clust_label <- c(1:ncol(ccc_unlist))
    colnames(ccc_unlist) <- clust_label


    clust_same <- (clust_cos > 0.95 & lower.tri(clust_cos))
    same <- which(clust_same, arr.ind=TRUE) # merge these columns
    if (length(same)>0){
      for (index in 1:nrow(same)){
        clust_label[same[index, 1]] <- clust_label[same[index, 2]]
      }
    }
    ccc_unlist <- data.frame(merge_cols(as.matrix(ccc_unlist),clust_label))
    summary <- c(summary,list(spectrum = ccc_unlist,
                              stats    = data.frame(table(clust_label))))

  }

  summary <- lapply(ccc_0,cosmergechain)

  dataframe <- data.frame(matrix(nrow=ncat,ncol=0))
  stats.dataframe <- data.frame(matrix(nrow=0,ncol=2))

  for(i in 1:nch){
    dataframe <- cbind(dataframe,summary[[i]]$spectrum)
    stats.dataframe <- rbind(stats.dataframe,summary[[i]]$stats)

  }


  dataframe <- dataframe[,stats.dataframe$Freq>1]##exclude spectrum that only extracted in 1 post sample in a chain
  stats.dataframe <- stats.dataframe[stats.dataframe$Freq>1,]

  clust_label <- 1:ncol(dataframe)
  colnames(dataframe) <- clust_label
  clust_label <- generate_label_high_cossim(clust_label = clust_label,
                                            matrix      = dataframe,
                                            cos.sim     = cos.merge)

  dataframe.merged <- data.frame(merge_cols(as.matrix(dataframe),clust_label))
  stats.dataframe.merged <- data.frame(merge_cols(t(stats.dataframe$Freq),clust_label))

  for(times in 1:5){
    clust_label <- c(1:ncol(dataframe.merged))
    colnames(dataframe.merged) <- clust_label

    names(stats.dataframe.merged) <- clust_label

    clust_cos <- cosCpp(as.matrix(dataframe.merged))

    clust_same <- (clust_cos > 0.90 & lower.tri(clust_cos))
    same <- which(clust_same, arr.ind=TRUE) # merge these columns
    if (length(same)>0){
      for (index in 1:nrow(same)){
        clust_label[same[index, 1]] <- clust_label[same[index, 2]]
      }
      dataframe.merged <- data.frame(merge_cols(as.matrix(dataframe.merged),clust_label))
      stats.dataframe.merged <- data.frame(merge_cols(as.matrix(stats.dataframe.merged),clust_label))
    }else{
      break
    }

  }
  dataframe.merged <- dataframe.merged[,order(colSums(dataframe.merged),decreasing=T)]
  stats.dataframe.merged <- stats.dataframe.merged[order(colSums(dataframe.merged),decreasing=T)]

  high.confident.spectrum <- dataframe.merged[,which(stats.dataframe.merged>=(0.9*nsamp))]
  high.confident.stats <- stats.dataframe.merged[which(stats.dataframe.merged>=(0.9*nsamp))]

  moderate.spectrum <- dataframe.merged[,intersect(which(stats.dataframe.merged>=(nsamp/nch)),which(stats.dataframe.merged<(0.9*nsamp)))]
  moderate.stats <- stats.dataframe.merged[intersect(which(stats.dataframe.merged>=(nsamp/nch)),which(stats.dataframe.merged<(0.9*nsamp)))]




  noise.spectrum <- dataframe.merged[,which(stats.dataframe.merged<(nsamp/nch))]
  noise.stats <- stats.dataframe.merged[which(stats.dataframe.merged<(nsamp/nch))]



  return(invisible(list(high.confident.spectrum = high.confident.spectrum,
                        high.confident.stats    = high.confident.stats,
                        moderate.spectrum = moderate.spectrum,
                        moderate.stats    = moderate.stats,
                        noise.spectrum = noise.spectrum,
                        noise.stats    = noise.stats,
                        ccc_0 = ccc_0,
                        cdc_0 = cdc_0,
                        multi.chains = x)))
}
