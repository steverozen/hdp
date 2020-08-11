#' Extract major components from the raw clusters
#'
#' If prior components included via \code{\link{hdp_prior_init}} are
#' preserved by \code{hdp_extract_components}, they are prefixed with "P".
#' Any new components in this case are prefixed with "N".
#'
#' @param x hdpSampleChain or hdpSampleMulti object
#' @param cos.merge Merge components with cosine similarity above this threshold (default 0.90)
#' @param min.sample Components must have significant exposure in at least this many samples (i.e. those DP nodes with data assigned) (default 1)
#' @param categ.CI A numeric between 0 and 1. Level of confidence interval to be calculated for each category.
#'                  Default is 0.95, but can be set to lower for extracting rare signatures
#' @param exposure.CI A numeric between 0 and 1. Level of confidence interval to be calculated
#'                    for a sample's exposure/observation of a raw cluster/proto-signature.
#'                    Default is 0.95, but can be set to lower for extracting rare signatures
#' @param diagnostic.folder If provided, details for hdp.0 is plotted
#'
#' @return A hdpSampleChain or hdpSampleMulti object updated with component information
#' @aliases hdp_merge_and_extract_components
#' @seealso \code{\link{hdp_posterior}}, \code{\link{hdp_multi_chain}},
#'  \code{\link{plot_comp_size}}, \code{\link{plot_comp_distn}},
#'  \code{\link{plot_dp_comp_exposure}}
#' @import clue
#' @export
# @examples
# hdp_extract_components(mut_example_multi)
hdp_merge_and_extract_components <- function(x,
                                             categ.CI    = 0.95,
                                             exposure.CI = 0.95,
                                             cos.merge   = 0.90,
                                             min.sample  = 1,
                                             diagnostic.folder=NULL){

  # input checks
  if (class(x)=="hdpSampleChain") {
    message('Extracting components on single chain.A hdpSampleMulti object is recommended, see ?hdp_multi_chain')
    is_multi <- FALSE
  } else if (class(x)=="hdpSampleMulti") {
    is_multi <- TRUE
  } else {
    stop("x must have class hdpSampleChain or hdpSampleMulti")
  }

  if (!validObject(x)) stop("x not a valid object")

  if (class(cos.merge) != "numeric" | cos.merge >=1 | cos.merge <= 0) {
    stop("cos.merge must be between 0 and 1")
  }
  if (min.sample %% 1 != 0 | min.sample < 1) {
    stop("min.sample must be a positive integer")
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


    #######################################

    ## At each sample point in the chain, merge clusters with cosine.similarity > 0.95
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

    ##every ccc and cdc has the same number of columns
    maxclust <- max(unlist((lapply(ccc_0, function(ccc){
      max(unlist(lapply(ccc,function(x)ncol(x))))
    } ))))

    for(i in 1:length(ccc_0)){
      ccc_temp <- ccc_0[[i]]
      cdc_temp <- cdc_0[[i]]

      for(j in 1:length(ccc_temp)){

        ccc_temp[[j]] <- cbind(ccc_temp[[j]],matrix(0, nrow=ncat, ncol=(maxclust-ncol(ccc_temp[[j]])+1)) )
        cdc_temp[[j]] <- cbind(cdc_temp[[j]],matrix(0, nrow=ndp, ncol=(maxclust-ncol(cdc_temp[[j]])+1)) )
      }
      ccc_0[[i]] <- ccc_temp
      cdc_0[[i]] <- cdc_temp

    }
    ###############################################################

    # if priors, #remove pseudo-counts from ccc_0

    ccc_raw_avg_per_ch <- lapply(ccc_0, function(matlist){ Reduce('+', matlist)/length(matlist) })

    mclust <- ncol(ccc_raw_avg_per_ch[[1]])

    rapch_unlist <- t(do.call(cbind, ccc_raw_avg_per_ch))
    rapch_gf <- rep(1:nch, each=mclust)
    rapch_ic <- rep(1:mclust, times=nch)
    ###re-indexing
    rapch_clust <- flexclust::kcca(rapch_unlist, k=rapch_ic,
                                   group=rapch_gf,
                                   family=flexclust::kccaFamily("kmedians",
                                                                groupFun="differentClusters"))

    rapch_label <- split(flexclust::clusters(rapch_clust), rapch_gf)

    ccc_1 <- Reduce('c', mapply(function(matlist, rank){
      lapply(matlist, function(mat){
        ans <- mat[,order(rank)]
        return(ans)
      })
    }, ccc_0, rapch_label, SIMPLIFY=FALSE))

    cdc_1 <- Reduce('c', mapply(function(matlist, rank){
      lapply(matlist, function(mat){
        ans <- mat[,order(rank)]
        return(ans)
      })
    }, cdc_0, rapch_label, SIMPLIFY=FALSE))

    remove(ccc_0, cdc_0, ccc_raw_avg_per_ch, rapch_unlist, rapch_gf, rapch_ic,
           rapch_clust, rapch_label, mclust)

  } else {

    maxclust <- max(numcluster(x))
    clust_label <- 1:maxclust

    ccc_1 <- lapply(clust_categ_counts(x), function(x){
      ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxclust-ncol(x)+1)))
      return(ans[, -ncol(ans)])
    })

    cdc_1 <- lapply(clust_dp_counts(x), function(x){
      ans <- cbind(x, matrix(0, nrow=ndp, ncol=(maxclust-ncol(x)+1)))
      return(ans[, -ncol(ans)])
    })

    # if priors, remove pseudo-counts from ccc_1
    if (is_prior){
      pseudodata <- sapply(dp(final_hdpState(x))[pseudo],
                           function(x) table(factor(x@datass, levels=1:ncat)))

      ccc_1 <- lapply(ccc_1, function(x) {
        x[,priorcc] <- x[,priorcc] - pseudodata
        return(x)
      })
    }

  }


  # Step (2)
  # Match up raw clusters (matrix columns) across posterior samples (columns not
  # guaranteed to keep same component through all samples)

  # K-centroids clustering of all raw clusters with cannot-link constraints
  # within each posterior sample, Manhattan distance and median centroid


  mclust <- ncol(ccc_1[[1]])

  if (mclust==1){
    ccc_label <- rep(1, length(ccc_1))

  } else{
    ccc_unlist <- t(do.call(cbind, ccc_1))

    # Change from NR code - the k medians clustering is
    # sensitive to the raw counts, so somethimes
    # clusters with the same profile are not identified
    # as such. So here we normalize by dividing the the
    # partial spectrum by the total number of mutations
    # in it.
    for(i in 1:nrow(ccc_unlist)){
      if(sum(ccc_unlist[i,])>0){
        ccc_unlist[i,] <- ccc_unlist[i,]/sum(ccc_unlist[i,])
      }
    }
    groupfactor <- rep(1:(nsamp), each=mclust)
    initial_clust <- rep(1:mclust, times=nsamp)

    ccc_clust <- flexclust::kcca(ccc_unlist, k=initial_clust,
                                 group=groupfactor,
                                 family=flexclust::kccaFamily("kmedians",

                                                              groupFun="differentClusters"))

    # want this plot to be as simple as possible
    # tmp <- matrix(flexclust::clusters(ccc_clust), byrow=T, ncol=mclust)
    # matplot(tmp, type='l', lty=1, main="kmedians")

    ccc_label <- split(flexclust::clusters(ccc_clust), groupfactor)

    remove(ccc_unlist, groupfactor, initial_clust, ccc_clust)

  }

  ccc_2 <- mapply(function(ccc, label) {
    colnames(ccc) <- label
    ccc[, order(as.numeric(colnames(ccc)))]
  },
  ccc_1, ccc_label, SIMPLIFY=FALSE)

  cdc_2 <- mapply(function(cdc, label) {
    colnames(cdc) <- label
    cdc[, order(as.numeric(colnames(cdc)))]
  },
  cdc_1, ccc_label, SIMPLIFY=FALSE)

  maxclust <- mclust
  clust_label <- 1:maxclust

  remove(ccc_1, cdc_1, ccc_label)

  # Step (3)
  # Merge the ccc columns with high cosine similarity.
  avgdistn <- matrix(0, nrow=ncat, ncol=maxclust)
  for (i in 1:maxclust){
    distns <- sapply(ccc_2, function(x) x[, i]/sum(x[, i]))
    avgdistn[, i] <- rowMeans(distns, na.rm=T)
  }
  clust_cos <- lsa::cosine(avgdistn)
  clust_same <- (clust_cos > cos.merge & lower.tri(clust_cos))
  same <- which(clust_same, arr.ind=TRUE) # merge these columns

  # update clust_label vector to reflect the merging of columns.
  if (length(same)>0){
    for (i in 1:nrow(same)){
      clust_label[same[i, 1]] <- clust_label[same[i, 2]]
    }
    #remove(i)
  }
  avgdistn_ccc3 <- merge_cols(avgdistn,clust_label)
  ccc_3 <- lapply(ccc_2, merge_cols, clust_label)
  cdc_3 <- lapply(cdc_2, merge_cols, clust_label)
  clust_label <- colnames(ccc_3[[1]])
  if (any(clust_label != colnames(cdc_3))) stop("problem in step 3!")

  ##more steps of doing merging. Some merged raw clusters have high cos sim
  for(iter.index in 1:10){

    avgdistn <- matrix(0, nrow=ncat, ncol=ncol(ccc_3[[1]]))
    for (i in 1:ncol(ccc_3[[1]])){
      distns <- sapply(ccc_3, function(x) x[, i]/sum(x[, i]))
      avgdistn[, i] <- rowMeans(distns, na.rm=T)
    }
    clust_cos <- lsa::cosine(avgdistn)
    clust_same <- (clust_cos > cos.merge & lower.tri(clust_cos))
    same <- which(clust_same, arr.ind=TRUE) # merge these columns
    if(length(same)==0){
      message("no more merging")
      break
    }else{
      message("extra merging")
      for (i in 1:nrow(same)){
        clust_label[same[i, 1]] <- clust_label[same[i, 2]]
      }
      #remove(i)
      avgdistn_ccc3 <- merge_cols(avgdistn,clust_label)
      ccc_3 <- lapply(ccc_3, merge_cols, clust_label)
      cdc_3 <- lapply(cdc_3, merge_cols, clust_label)
      clust_label <- colnames(ccc_3[[1]])
      if (any(clust_label != colnames(cdc_3))) stop("problem in step 3!")
    }
    # update clust_label vector to reflect the merging of columns.

  }


  remove(avgdistn, distns, clust_cos, clust_same, same, ccc_2, cdc_2)


  # Step (4)
  # Assign components with no *significantly* non-zero data categories
  # to component '0'
  clust_hdp0_ccc4 <- data.frame(matrix(ncol=0,nrow=ncat))
  use_clust <- c()
  for (ii in 1:ncol(ccc_3[[1]])) {
    compii <- sapply(ccc_3, function(x) x[,ii])
    lowerb <- apply(compii, 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        NaN
      } else {
        round(coda::HPDinterval(samp, categ.CI)[1], 3)

      }
    })
    if(any(lowerb>0)){
      use_clust <- c(use_clust, colnames(ccc_3[[1]])[ii])
    }else{
      ##I want to summarize how many mutations of which cluster were moved into hdp.0 -- Mo
      clust_hdp0_ccc4 <- cbind(clust_hdp0_ccc4,rowSums(compii))

      colnames(clust_hdp0_ccc4)[ncol(clust_hdp0_ccc4)] <- paste("ccc_3_",colnames(ccc_3[[1]])[ii],sep="")
    }
  }


  # update clust_label vector
  clust_label[which(!clust_label %in% use_clust)] <- '0'

  ccc_4 <- lapply(ccc_3, merge_cols, clust_label)

  cdc_4 <- lapply(cdc_3, merge_cols, clust_label)

  # Change from NR code: back pointers to the constuents
  # ccc_4 element.
  avgdistn_ccc4 <- matrix(0, nrow=ncat, ncol=ncol(ccc_4[[1]]))

  for (i in 1:ncol(ccc_4[[1]])){
    distns <- sapply(ccc_4, function(x) x[, i])
    avgdistn_ccc4[, i] <- rowSums(distns, na.rm=T)
  }
  colnames(avgdistn_ccc4) <- colnames(cdc_4[[1]])
  clust_label <- colnames(cdc_4[[1]])

  if (any(clust_label != colnames(cdc_4))) stop("problem in step 4!")


  # if there was no component zero added, add an empty one now
  if (!"0" %in% clust_label) {
    ccc_4 <- lapply(ccc_4, function(x){
      ans <- cbind(0, x)
      colnames(ans) <- c(0, colnames(x))
      return(ans)
    })


    cdc_4 <- lapply(cdc_4, function(x){
      ans <- cbind(0, x)
      colnames(ans) <- c(0, colnames(x))
      return(ans)
    })

  }

  ##Added by Mo: create diagnostic plot to trace back hdp.0

  if(!is.null(diagnostic.folder)){
    if (dir.exists(diagnostic.folder)) {
      message(diagnostic.folder, " already exits")
    } else {
      dir.create(diagnostic.folder, recursive = T)
    }

    ##plot1

    row.names(avgdistn_ccc4) <- ICAMS::catalog.row.order$SBS96
    avgdistn_ccc4_catalog <- ICAMS::as.catalog(avgdistn_ccc4,catalog.type = "counts")
    ICAMS::PlotCatalogToPdf(avgdistn_ccc4_catalog,
                            file.path(diagnostic.folder, "aggregated.spectrum.after.step4.pdf"))


    if(ncol(clust_hdp0_ccc4)>0&&!is.null(ncol(clust_hdp0_ccc4))){

      ##plot2

      row.names(clust_hdp0_ccc4) <- ICAMS::catalog.row.order$SBS96
      clust_hdp0_ccc4_catalog <- ICAMS::as.catalog(clust_hdp0_ccc4,catalog.type = "counts")
      ICAMS::PlotCatalogToPdf(clust_hdp0_ccc4_catalog,
                              file.path(diagnostic.folder, "aggregated.spectrum.moved.to.hdp0.in.step4.pdf"))

      ##plot3 each cluster in clust_hdp0_ccc4 has a new folder


      diagnostic_in_extraction(clust_hdp0_ccc = clust_hdp0_ccc4,
                               nsamp          = nsamp,
                               ncat           = ncat,
                               nch            = nch,
                               ccc            = ccc_3,
                               cdc            = cdc_3,
                               diagnostic.folder = diagnostic.folder)


    }
  }


  remove(compii, ccc_3, cdc_3, ii, lowerb, use_clust)

  # Step (5)
  # Assign components with < min.sample *significantly* non-zero sample exposure
  # to component '0' (disregarding DP nodes with no data items (parent nodes))
  use_clust <- c()
  disregard <- if(is_prior) union(which(numdata==0), pseudo) else which(numdata==0)
  clust_hdp0_ccc5 <- data.frame(matrix(ncol=0,nrow=ncat))

  for (ii in 1:ncol(cdc_4[[1]])) {
    compii <- sapply(cdc_4, function(x) x[,ii])
    lowerb <- apply(compii[-disregard,], 1, function(y) {

      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        NaN
      } else {
        round(coda::HPDinterval(samp, exposure.CI)[1], 3)

      }
    })
    if(sum(lowerb>0)>=min.sample){
      use_clust <- c(use_clust, colnames(cdc_4[[1]])[ii])
    }else{
      ##I want to summarize how many mutations of which cluster were moved into hdp.0 -- Mo
      ccc_compii <- sapply(ccc_4, function(x) x[,ii])

      clust_hdp0_ccc5 <- cbind(clust_hdp0_ccc5,rowSums(ccc_compii))

      colnames(clust_hdp0_ccc5)[ncol(clust_hdp0_ccc5)] <- paste("ccc_4_",colnames(ccc_4[[1]])[ii],sep="")
    }
  }


  # update clust_label vector
  clust_label[which(!clust_label %in% use_clust)] <- 0
  ccc_5 <- lapply(ccc_4, merge_cols, clust_label)
  cdc_5 <- lapply(cdc_4, merge_cols, clust_label)
  clust_label <- colnames(ccc_5[[1]])
  if (any(clust_label != colnames(cdc_5))) stop("problem in step 5!")


  ##Added by Mo: create diagnostic plot to trace back hdp.0

  if(!is.null(diagnostic.folder)){
    if (dir.exists(diagnostic.folder)) {
      message(diagnostic.folder, " already exits")
    } else {
      dir.create(diagnostic.folder, recursive = T)
    }

    ##plot1

    row.names(avgdistn_ccc4) <- ICAMS::catalog.row.order$SBS96
    avgdistn_ccc4_catalog <- ICAMS::as.catalog(avgdistn_ccc4,catalog.type = "counts")
    ICAMS::PlotCatalogToPdf(avgdistn_ccc4_catalog,file.path(diagnostic.folder, "aggregated.spectrum.after.step5.pdf"))

    if(ncol(clust_hdp0_ccc5)>0&&!is.null(ncol(clust_hdp0_ccc5))){

      ##plot2

      row.names(clust_hdp0_ccc5) <- ICAMS::catalog.row.order$SBS96
      clust_hdp0_ccc5_catalog <- ICAMS::as.catalog(clust_hdp0_ccc5,catalog.type = "counts")
      ICAMS::PlotCatalogToPdf(clust_hdp0_ccc5_catalog,
                              file.path(diagnostic.folder, "aggregated.spectrum.moved.to.hdp0.in.step5.pdf"))

      ##plot3 each cluster in clust_hdp0_ccc4 has a new folder


      diagnostic_in_extraction(clust_hdp0_ccc = clust_hdp0_ccc5,
                               nsamp          = nsamp,
                               ncat           = ncat,
                               nch            = nch,
                               ccc            = ccc_4,
                               cdc            = cdc_4,
                               diagnostic.folder = diagnostic.folder)
    }
  }




  remove(compii, ccc_4, cdc_4, ii, lowerb, use_clust, disregard)

  # Step (6)
  # Rename overall component, order by number of data items (on average)
  # 0th component still goes first

  avg_ndi <- rowMeans(sapply(ccc_5, colSums))
  colorder <- c(1, setdiff(order(avg_ndi, decreasing=T), 1))


  # If priors,
  # update clust_label to reflect match (down to 0.9) with prior components
  if (is_prior) {

    nco <- length(clust_label)

    # average distribution over data categ for each component
    avgdistn <- sapply(2:nco, function(i){
      distns <- sapply(ccc_5, function(x) x[, i]/sum(x[, i]))
      ans <- rowMeans(distns, na.rm=T)
      return(ans)
    })

    # compare against original pseudodata distn
    match2_pseudo <- apply(avgdistn, 2, function(x) lsa::cosine(x, pseudodata))
    rownames(match2_pseudo) <- priorcc
    colnames(match2_pseudo) <- clust_label[-1]

    to_match <- TRUE
    while(to_match){
      if(!any(match2_pseudo>0.9)) {
        to_match <- FALSE
        break
      }

      best_match <- which(match2_pseudo==max(match2_pseudo), arr.ind = TRUE)
      old <- colnames(match2_pseudo)[best_match[2]]
      new <- paste0("P", rownames(match2_pseudo)[best_match[1]])
      clust_label[which(clust_label == old)] <- new

      match2_pseudo <- match2_pseudo[-best_match[1], -best_match[2], drop=FALSE]

    }
    suppressWarnings(rm(to_match, match2_pseudo, avgdistn, best_match, old, new, nco))

    ccc_5 <- lapply(ccc_5, function(x) {
      colnames(x) <- clust_label
      return(x)
    })

    cdc_5 <- lapply(cdc_5, function(x) {
      colnames(x) <- clust_label
      return(x)
    })

  }



  ccc_6 <- lapply(ccc_5, function(x) {
    x <- x[, colorder]
    if (is_prior){
      update <- setdiff(which(!grepl("P", colnames(x))), 1)
      if (length(update)>0){
        colnames(x)[update] <- paste0("N", 1:length(update))
      }
    } else {
      colnames(x) <- 0:(ncol(x)-1)
    }
    return(x)
  })

  cdc_6 <- lapply(cdc_5, function(x) {
    x <- x[, colorder]
    if (is_prior){
      update <- setdiff(which(!grepl("P", colnames(x))), 1)
      if (length(update)>0){
        colnames(x)[update] <- paste0("N", 1:length(update))
      }
    } else {
      colnames(x) <- 0:(ncol(x)-1)
    }
    return(x)
  })

  # number of components
  ncomp <- length(colorder)

  remove(ccc_5, cdc_5, avg_ndi, colorder)


  # Step (7)
  # Convert ccc into list of length ncomp, with matrices nsamp*ncat

  ccc_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncat)), ncomp)
  for (i in 1:ncomp){
    ccc_ans[[i]] <- t(sapply(ccc_6, function(x) x[, i]))
  }
  names(ccc_ans) <- colnames(ccc_6[[1]])

  # Convert cdc into list of length ndp, with matrices nsamp*ncomp
  cdc_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncomp)), ndp)
  for (i in 1:ndp){
    cdc_ans[[i]] <- t(sapply(cdc_6, function(x) x[i, ]))
  }

  remove(ccc_6, cdc_6)

  # Step (8)
  # Calculate mean and 95% credibility interval for each component's
  # categorical data distribution
  ccc_norm <- lapply(ccc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  ccc_mean <- t(sapply(ccc_norm, colMeans, na.rm=TRUE))

  ccc_credint <- lapply(ccc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, 0.95), 4)
      }
    })
  })

  # Step (8)
  # Calculate mean and 95% credibility interval for each DP's
  # distribution over components (counts)
  cdc_norm <- lapply(cdc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  cdc_mean <- t(sapply(cdc_norm, colMeans, na.rm=TRUE))

  cdc_credint <- lapply(cdc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, 0.95), 4)
      }
    })
  })

  # add extracted components into x hdpSampleChain slots
  x@numcomp <- as.integer(ncomp - 1)

  # proportion of data explained by extracted components?
  avcount <- colMeans(sapply(ccc_ans, rowSums, na.rm=TRUE), na.rm=TRUE)
  x@prop.ex <- round(1-avcount[1]/sum(avcount), 3)

  x@comp_cos_merge <- cos.merge


  x@comp_categ_counts <- ccc_ans
  x@comp_dp_counts <- lapply(cdc_ans, as, "dgCMatrix")


  x@comp_categ_distn <- list(mean                                        = ccc_mean,
                             cred.int                                    = ccc_credint,
                             aggregated_raw_clusters_after_cos_merge     = avgdistn_ccc3,
                             aggregated_raw_clusters_after_nonzero_categ = avgdistn_ccc4,
                             clust_hdp0_ccc4                             = clust_hdp0_ccc4,
                             clust_hdp0_ccc5                             = clust_hdp0_ccc5)

  x@comp_dp_distn <- list(mean=cdc_mean,
                          cred.int=cdc_credint)

  # check validity and return
  if (!validObject(x)) warning("Not a valid hdpSampleChain/Multi object.")
  return(x)
}
