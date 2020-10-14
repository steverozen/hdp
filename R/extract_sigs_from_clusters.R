#' This function extract signatures from raw clusters discovered by hdp.
#' This returns clusters with the number of posterior samples that they belong
#'
#' ccc_0 and cdc_0 are for diagnostic checking and exposure summary
#'
#' @param x \code{hdpSampleChain} or \code{hdpSampleMulti} object
#' @param cos.merge Merge components with cosine similarity above this threshold (default 0.90).
#' @param hc.cutoff the height to cut hierarchcial clustering tree
#' @return A hdpSampleChain or hdpSampleMulti object updated with component information
#' @aliases extract_sigs_from_clusters
#' @seealso \code{\link{hdp_posterior}}, \code{\link{hdp_multi_chain}},
#'  \code{\link{plot_comp_size}}, \code{\link{plot_comp_distn}},
#'  \code{\link{plot_dp_comp_exposure}}
#' @importFrom stats cutree aggregate

#' @export


extract_sigs_from_clusters <-  function(x,
                                        cos.merge = 0.90,
                                        hc.cutoff = 0.12
                                      ){
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
    clust_same <- (clust_cos > 0.99 & lower.tri(clust_cos))
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


  for(i in 1:nch){

    test <- mapply(first.merge,ccc_0[[i]],cdc_0[[i]],SIMPLIFY = F)

    for(j in 1:length(ccc_0[[i]])){
      ccc_0[[i]][[j]] <- test[[j]]$ccc
      cdc_0[[i]][[j]] <- test[[j]]$cdc
    }



  }
  ########################################################################
  #######merge clusters with high cos.sim in one posterior chain##########
  ########################################################################


  cosmergechain <- function(ccc,cdc){
    clust_cos <- cosCpp(as.matrix(ccc))
    clust_label <- c(1:ncol(ccc))

    colnames(ccc) <- colnames(cdc) <- clust_label

    clust_same <- (clust_cos > 0.99 & lower.tri(clust_cos))
    same <- which(clust_same, arr.ind=TRUE) # merge these columns
    if (length(same)>0){
      for (index in 1:nrow(same)){
        clust_label[same[index, 1]] <- clust_label[same[index, 2]]
      }
    }


    ccc_unlist <- merge_cols(as.matrix(ccc),clust_label)
    cdc_unlist <- merge_cols(as.matrix(cdc),clust_label)

    stats <- data.frame(table(clust_label))

    new.spectrum <- ccc_unlist[,which(stats$Freq>1)]
    new.cdc <- cdc_unlist[,which(stats$Freq>1)]
    new.stats <- stats[which(stats$Freq>1),]
    noise.spectrum <- rowSums(ccc_unlist[,which(stats$Freq==1)])
    noise.cdc <- rowSums(cdc_unlist[,which(stats$Freq==1)])
    noise.stats <- sum(stats$Freq==1)
    return(list= list(spectrum     = new.spectrum,
                      spectrum_cdc = new.cdc,
                      spectrum_stats = new.stats,
                      noise.spectrum = noise.spectrum,
                      noise.cdc    = noise.cdc,
                      noise.stats  = noise.stats
                      ))

  }

  for(i in 1:nch){
    ccc_0[[i]] <-  do.call(cbind,ccc_0[[i]])
    cdc_0[[i]] <- do.call(cbind,cdc_0[[i]])
  }
  summary <- mapply(cosmergechain,ccc_0,cdc_0,SIMPLIFY = F)

  dataframe <- all.noise.spectrum <- data.frame(matrix(nrow=ncat,ncol=0))
  stats.dataframe <- data.frame(matrix(nrow=0,ncol=2))
  dp.dataframe <- each.chain.noise.cdc <- data.frame(matrix(nrow=ndp,ncol=0))

  for(i in 1:nch){
    dataframe <- cbind(dataframe,summary[[i]]$spectrum)
    dp.dataframe <- cbind(dp.dataframe,summary[[i]]$spectrum_cdc)
    stats.dataframe <- rbind(stats.dataframe,summary[[i]]$spectrum_stats)
    each.chain.noise.spectrum <- cbind(all.noise.spectrum,summary[[i]]$noise.spectrum)
    each.chain.noise.cdc <- cbind(each.chain.noise.cdc,summary[[i]]$noise.cdc)
  }


  #dataframe <- dataframe[,stats.dataframe$Freq>1]##exclude spectrum that only extracted in 1 post sample in a chain, too many noisy clusters affect the final extraction
  #stats.dataframe <- stats.dataframe[stats.dataframe$Freq>1,]

  dataframe.normed <- apply(dataframe,2,function(x)x/sum(x))
  cosine.dist.df <- parallelDist::parallelDist(t(dataframe.normed),method = "cosine")
  cosine.dist.hctree <- stats::hclust(cosine.dist.df,method = "average")

  ####decide best cutoff##########################################
  #####cut hc tree until no clusters have cos.sim > cos.cutoff####

  clusters <- cutree(cosine.dist.hctree,  h=hc.cutoff) #make sure each cluster is clean
  spectrum.df <- merge_cols(as.matrix(dataframe),clusters)
  spectrum.stats <- aggregate(stats.dataframe[,2],by=list(clusters),sum)
  spectrum.cdc <- merge_cols(as.matrix(dp.dataframe),clusters)


  for(iter.index in 1:15){

    clust_cos <- cosCpp(as.matrix(spectrum.df))
    clust_label <- c(1:ncol(spectrum.df))
    colnames(spectrum.df) <- c(1:ncol(spectrum.df))
    colnames(spectrum.cdc) <- c(1:ncol(spectrum.df))
    clust_same <- (clust_cos > cos.merge & lower.tri(clust_cos))
    same <- which(clust_same, arr.ind=TRUE) # merge these columns
    if(length(same)==0){
      #message("no more merging")
      break
    }else{
      #message("extra merging")
      for (i in 1:nrow(same)){
        clust_label[same[i, 1]] <- clust_label[same[i, 2]]
      }
      #remove(i)
      spectrum.df <- merge_cols(as.matrix(spectrum.df),clust_label)
      spectrum.stats <- aggregate(spectrum.stats[,2],by=list(clust_label),sum)
      spectrum.cdc <- merge_cols(as.matrix(spectrum.cdc),clust_label)
    }
    # update clust_label vector to reflect the merging of columns.


  }

  return(invisible(list(clustered.spectrum = spectrum.df,
                        stats.post.samples = spectrum.stats,
                        spectrum.cdc = spectrum.cdc,
                        each.chain.noise.cdc = each.chain.noise.cdc,  #noise from each chain

                        each.chain.noise.spectrum = each.chain.noise.spectrum,  #noise from each chain

                        multi.chains = x,
                        nsamp = nsamp)))
}
