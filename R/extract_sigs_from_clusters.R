#' This function combines part of NR's code and new way to look at sigantures from raw clusters
#' This is a testing function that maybe or not work
#'
#' @param x hdpSampleChain or hdpSampleMulti object
#' @param cos.merge Merge components with cosine similarity above this threshold (default 0.90)
#' @param min.sample Components must have significant exposure in at least this many samples (i.e. those DP nodes with data assigned) (default 1)
#' @param prop.samp The proportion of posterior samples with a certain cluster to be determined as                     a signature. e.g., if prop.samp = 0.5, a cluster will be determined as a                           signature if it appears in half of the posterior samples
#'
#' @return A hdpSampleChain or hdpSampleMulti object updated with component information
#' @aliases extract_sigs_from_clusters
#' @seealso \code{\link{hdp_posterior}}, \code{\link{hdp_multi_chain}},
#'  \code{\link{plot_comp_size}}, \code{\link{plot_comp_distn}},
#'  \code{\link{plot_dp_comp_exposure}}
#' @importFrom stats sd
#' @import clue
#' @export
# @examples
# hdp_extract_components(mut_example_multi)
extract_sigs_from_clusters <- function(x,
                                       cos.merge,
                                       prop.samp,
                                       min.sample){
  # input checks
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

      len <- length(ccc_0[[i]])
      ccc_temp <- ccc_0[[i]]
      cdc_temp <- cdc_0[[i]]

      for(j in 1:length(ccc_temp)){

        ccc_temp[[j]] <- cbind(ccc_temp[[j]],matrix(0, nrow=ncat, ncol=(maxclust+1-ncol(ccc_temp[[j]]))) )
        cdc_temp[[j]] <- cbind(cdc_temp[[j]],matrix(0, nrow=ndp, ncol=(maxclust+1-ncol(cdc_temp[[j]]))) )
      }

      ccc_0[[i]] <- ccc_temp
      cdc_0[[i]] <- cdc_temp
      ####a first step of re-indexing?###



    }
    ###############################################################

    # if priors, #remove pseudo-counts from ccc_0

    ccc_raw_avg_per_ch <- lapply(ccc_0, function(matlist){ Reduce('+', matlist)/length(matlist) })

    mclust <- ncol(ccc_raw_avg_per_ch[[1]])

    rapch_unlist <- t(do.call(cbind, ccc_raw_avg_per_ch))

    for(index in 1:nrow(rapch_unlist)){
      if(sum(rapch_unlist[index,])>0){
        rapch_unlist[index,] <- rapch_unlist[index,]/sum(rapch_unlist[index,])
      }
    }

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
  remove(ccc_1, cdc_1, ccc_label)


  results.compiled <- {}
  disregard <- if(is_prior) union(which(numdata==0), pseudo) else which(numdata==0)

  #  parallel.run.identify.individual.cluster <- function(x){
  #    retval <- identify_individual_cluster(x,
  #                                          ccc_2,
  #                                          cdc_2,
  #                                          disregard)
  #return(retval)
  # }

  #  results.compiled <- parallel::mclapply(
  # Must choose a different seed for each of the chains
  #    X = 1:ncol(ccc_2[[1]]),
  #   FUN = parallel.run.identify.individual.cluster,
  #  mc.cores = ncol(ccc_2[[1]]))

  ###does mclapply makes this faster?


  for(i in 1:ncol(ccc_2[[1]])){

    retval <- identify_individual_cluster(i,
                                          ccc_2,
                                          cdc_2,
                                          disregard)


    results.compiled <- append(results.compiled, list(retval))
    names(results.compiled)[i] <-  paste0("potential.hdp.",i)

  }##time consuming

  results.compiled <- results.compiled[which(unlist(lapply(results.compiled,length))>0)]
  save(results.compiled,file="results.compiled.Rdata")
  ##Step4 selected signatures based on how many posterior samples having them

  selected.sig <-   nonselected.sig <-   selected.cdc <-   nonselected.cdc <-
    selected.ccc <-   nonselected.ccc <- selected.matched.samples <- nonselected.matched.samples<-{}



  for(i in 1:length(results.compiled)){
    results.temp <- results.compiled[i][[1]]
    for(j in 1:length(results.temp)){

      if(results.temp[[j]]$matched.samp>=nsamp*prop.samp){
        selected.sig <- cbind(selected.sig,results.temp[[j]]$counts.spec)
        selected.cdc <- append(selected.cdc,
                               list(results.temp[[j]]$exp.df))
        selected.ccc <- append(selected.ccc,
                               list(results.temp[[j]]$ccc.agg))
        selected.matched.samples <- c(selected.matched.samples,results.temp[[j]]$matched.samp)


      }else{
        nonselected.sig <- cbind(nonselected.sig,results.temp[[j]]$counts.spec)
        nonselected.cdc <- append(nonselected.cdc,list(results.temp[[j]]$exp.df))
        nonselected.ccc <- append(nonselected.ccc,
                                  list(results.temp[[j]]$ccc.agg))
        nonselected.matched.samples <- c(nonselected.matched.samples,results.temp[[j]]$matched.samp)
      }

    }


  }


  ##remove duplicated sigs
  clust_label <- 1:ncol(selected.sig)
  colnames(selected.sig) <- clust_label

  clust_label <- generate_label_high_cossim(clust_label = clust_label,
                                            matrix      = selected.sig,
                                            cos.sim     = 0.95)

  selected.sig <- selected.sig[,unique(clust_label)]
  selected.cdc <- selected.cdc[unique(clust_label)]
  selected.ccc <- selected.ccc[unique(clust_label)]
  selected.matched.samples <- selected.matched.samples[unique(clust_label)]


  clust_label <- 1:ncol(nonselected.sig)
  colnames(nonselected.sig) <- clust_label

  clust_label <- generate_label_high_cossim(clust_label = clust_label,
                                            matrix      = nonselected.sig,
                                            cos.sim     = 0.95)


  nonselected.sig <- nonselected.sig[,unique(clust_label)]
  nonselected.cdc <- nonselected.cdc[unique(clust_label)]
  nonselected.ccc <- nonselected.ccc[unique(clust_label)]
  nonselected.matched.samples <- nonselected.matched.samples[unique(clust_label)]



  ##merge similar sigs


  clust_label <- 1:ncol(selected.sig)
  colnames(selected.sig) <- clust_label
  clust_label <- generate_label_high_cossim(clust_label = clust_label,
                                            matrix      = selected.sig,
                                            cos.sim     = cos.merge)




  for(i in 1:length(clust_label)){
    if(length(which(clust_label==clust_label[i]))>1){

      selected.cdc[[i]] <- do.call(cbind,selected.cdc[index])
      selected.ccc[[i]] <- do.call(cbind,selected.ccc[index])
    }

  }

  selected.cdc <- selected.cdc[unique(clust_label)]
  selected.ccc <- selected.ccc[unique(clust_label)]

  selected.sig <- data.frame(merge_cols(as.matrix(selected.sig),clust_label))
  selected.matched.samples <- merge_cols(t(selected.matched.samples),clust_label)





  clust_label <- 1:ncol(nonselected.sig)
  colnames(nonselected.sig) <- clust_label
  clust_label <- generate_label_high_cossim(clust_label = clust_label,
                                            matrix      = nonselected.sig,
                                            cos.sim     = cos.merge)


  for(i in 1:length(clust_label)){
    if(length(which(clust_label==clust_label[i]))>1){
      index <- which(clust_label==clust_label[i])

      nonselected.cdc[[i]] <- do.call(cbind,nonselected.cdc[index])
      nonselected.ccc[[i]] <- do.call(cbind,nonselected.ccc[index])

    }

  }

  nonselected.cdc <- nonselected.cdc[unique(clust_label)]
  nonselected.ccc <- nonselected.ccc[unique(clust_label)]
  nonselected.sig <- nonselected.sig[,unique(clust_label)]
  nonselected.matched.samples <- merge_cols(t(nonselected.matched.samples),clust_label)


  for(i in 1:ncol(selected.sig)){
    compii <- selected.cdc[[i]]
    lowerb <- apply(compii, 1, function(y) {

      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        NaN
      } else {
        round(coda::HPDinterval(samp, 0.95)[1], 3)

      }
    })
    if(sum(lowerb>0)<min.sample){

      ##I want to summarize how many mutations of which cluster were moved into hdp.0 -- Mo
      nonselected.sig <- cbind(nonselected.sig,selected.sig[,i])

      nonselected.cdc <- c(nonselected.cdc,selected.cdc[[i]])
      nonselected.ccc <- c(nonselected.ccc,selected.ccc[[i]])
      nonselected.matched.samples <- c(nonselected.matched.samples,selected.matched.samples[[i]])


      selected.sig <- selected.sig[,-i]
      selected.cdc <- selected.cdc[-i]
      selected.ccc <- selected.ccc[,-i]
      selected.matched.samples <- selected.matched.samples[-i]


    }

  }

  ###export individual clusters that with samp.prop less than required
  ###useful for a future selected
  noise.individual.clusters <- nonselected.sig
  colnames(noise.individual.clusters) <- paste0("noise.cluster.",1:ncol(noise.individual.clusters))
  names(nonselected.matched.samples) <- colnames(noise.individual.clusters)
  order1 <- order(nonselected.matched.samples,decreasing = T)
  noise.individual.clusters <- noise.individual.clusters[,order1]
  nonselected.matched.samples <- nonselected.matched.samples[,order1]
  nonselected.cdc <- do.call(cbind,nonselected.cdc)
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################

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


  nonselected.sig.sum <- rowSums(nonselected.sig)##generate hdp.0 by sum of all nonselected.sigs
  selected.sig <- cbind(nonselected.sig.sum,selected.sig)
  colorder <- c(1, setdiff(order(colSums(selected.sig), decreasing=T), 1))
  selected.sig <- selected.sig[,colorder]

  # number of components
  ncomp <- length(colorder)

  selected.sig <- apply(selected.sig,2,function(x){x/sum(x)})

  for(i in 1:length(selected.cdc)){
    selected.cdc[[i]] <- do.call(cbind,selected.cdc[[i]])
  }

  selected.exp.means <- lapply(selected.cdc,rowMeans)
  nonselected.exp.means <- rowMeans(nonselected.cdc)

  selected.exp <- do.call(cbind, selected.exp.means)

  comp_all_exp <- cbind(nonselected.exp.means,selected.exp)

  comp_all_exp <- comp_all_exp[,colorder]

  comp_all_exp <- t(apply(comp_all_exp,1,function(x){x/sum(x)}))

  #######

  aggregated.nonselected.ccc <- do.call(cbind, nonselected.ccc)

  all.ccc <- c(list(aggregated.nonselected.ccc),selected.ccc)

  max.nsamp <- max(sapply(all.ccc, function(x) ncol(x)))

  aggregated.nonselected.cdc <- do.call(cbind, nonselected.cdc)

  all.cdc <- c(list(aggregated.nonselected.cdc),selected.cdc)


  cdc_ans <- rep(list(matrix(0, nrow=max.nsamp, ncol=ncomp)), (ndp-length(disregard)))


  for (i in 1:(ndp-length(disregard))){

    temp <-t(lapply(all.cdc, function(x) x[i, ]))

    for(len in 1:length(temp)){
      cdc_ans[[i]][,len] <- c(as.numeric(unlist(temp[[len]])),rep(0,max.nsamp-length(temp[[len]])))
    }

  }

  ccc_norm <- lapply(all.ccc, function(x)  sweep(x, 2, colSums(x), FUN="/"))

  ccc_mean <- sapply(ccc_norm, rowMeans, na.rm=TRUE)
  colnames(ccc_mean) <- colorder

  ccc_credint <- lapply(ccc_norm, function(x) {
    apply(x, 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, 0.95), 4)
      }
    })
  })
  names(ccc_credint) <- colorder

  # Step (8)
  # Calculate mean and 95% credibility interval for each DP's
  # distribution over components (counts)
  cdc_norm <- lapply(all.cdc, function(x) sweep(x, 2, colSums(x), FUN="/"))

  cdc_mean <- sapply(cdc_norm, rowMeans, na.rm=TRUE)

  cdc_credint <- lapply(cdc_norm, function(x) {
    apply(x, 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, 0.95), 4)
      }
    })
  })

  all.ccc <- lapply(all.ccc,t)  ##transform to fit for diagnostic plotting



  # add extracted components into x hdpSampleChain slots
  x@numcomp <- as.integer(ncomp - 1)

  #avcount <- mean(sapply(all.ccc, rowSums, na.rm=TRUE), na.rm=TRUE)

  x@prop.ex <- round(1-colSums(comp_all_exp)[1]/sum(colSums(comp_all_exp)), 3)

  x@comp_cos_merge <- cos.merge


  x@comp_categ_counts <- all.ccc

  x@comp_dp_counts <- lapply(cdc_ans, as, "dgCMatrix")



  x@comp_cos_merge <- cos.merge

  x@comp_categ_distn <- list(mean             = t(ccc_mean),
                             cred.int         = ccc_credint,
                             comp_categ_distn = selected.sig,
                             noise.individual.clusters = noise.individual.clusters,
                             nonselected.matched.samples = nonselected.matched.samples)

  x@comp_dp_distn <- list(mean=comp_all_exp,
                          cred.int=cdc_credint)

  # check validity and return
  #if (!validObject(x)) warning("Not a valid hdpSampleChain/Multi object.")
  return(x)
}




