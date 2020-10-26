#'find the ccc and cdc that matched to a spectrum in ccc_0 and cdc_0
#'this function is to summarize the credint and mean of cccs and cdcs
#'and for further diagnostic plotting
#'@param spectrum signature spectrum for compare
#'@param ccc_0 a list object contains clust_categ_counts matrix from hdp
#'@param cdc_0 a list object contains clust_dp_counts matrix from hdp
#'@param cos.merge cosine similarity cutoff
#'
#'@export
extract_ccc_cdc_from_hdp <- function(spectrum,
                                     ccc_0,
                                     cdc_0,
                                     cos.merge = 0.90){

  spectrum.ccc <- data.frame(matrix(nrow=nrow(ccc_0[[1]][[1]]),ncol=0))
  summary.chain.info <- data.frame(matrix(ncol=3,nrow=0))



  for(chain in 1:length(ccc_0)){

    temp.chain <- data.frame(matrix(nrow=length(ccc_0[[chain]]),ncol=3))
    temp.chain[,1] <- chain
    temp.chain[,2] <-  temp.chain[,4] <- 0
    temp.chain[,3] <- c(1:nrow(temp.chain))

    for(sample in 1:length(ccc_0[[chain]])){

      ccc_0_temp <- ccc_0[[chain]][[sample]]

      cos.sims <- apply(ccc_0_temp,2,function(x){lsa::cosine(x,spectrum)})
      if(sum(cos.sims>cos.merge)>0){
        spectrum.ccc <- cbind(spectrum.ccc,ccc_0_temp[,cos.sims>cos.merge])
        temp.chain[sample,2] <- chain
        temp.chain[sample,4] <- sum(colSums(ccc_0_temp[,cos.sims>cos.merge,drop=F]))

      }

    }
    summary.chain.info <- rbind(summary.chain.info,temp.chain)


  }

  ccc_norm <- apply(spectrum.ccc, 2,function(x) x/sum(x, na.rm=TRUE))

  ccc_mean <- rowMeans(ccc_norm)

  ccc_credint <- apply(ccc_norm, 1, function(y) {
    samp <- coda::as.mcmc(y)
    if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
      c(NaN, NaN)
    } else {
      round(coda::HPDinterval(samp, 0.95), 4)
    }
  })


  return(invisible(list(spectrum.ccc = spectrum.ccc,
                        ccc_mean = ccc_mean,
                        ccc_credint = ccc_credint,
                       summary.chain.info=summary.chain.info)))
}

#' Plot signatures and their 95\% credible intervals
#'
#' @param retval an object return from \code{\link{extract_ccc_cdc_from_hdp}}.
#' @param col Either a single colour for all data categories, or a vector of
#'  colours for each group (in the same order as the levels of the grouping factor)
#' @param cred_int Logical - should 95\% credibility intervals be plotted? (default TRUE)
#' @param weights (Optional) Weights over the data categories to adjust their
#'  relative contribution (multiplicative)

#' @param group_label_height Multiplicative factor from top of plot for group label placement
#' @param cat_names names displayed on x-axis, e.g. SBS96 mutation classes
#' @param cex.cat Expansion factor for the (optional) cat_names
#'
#' @export
#'
plot_component_with_credint <-
  function(retval,cat_names=NULL,
           col="grey70",
           cred_int=TRUE,
           weights=NULL,
           group_label_height=1.05, cex.cat=0.7){

    # input checks
    ccc_mean_df <- do.call(cbind,lapply(retval,function(x)x[["ccc_mean"]]))
    ccc_credint <- lapply(retval,function(x)x[["ccc_credint"]])

    ncat <- nrow(ccc_mean_df)

    comp_distn <- ccc_mean_df

    if(class(cred_int) != "logical") stop("cred_int must be TRUE or FALSE")
    if(!class(weights) %in% c("numeric", "NULL") |
       !length(weights) %in% c(ncat, 0)) {
      stop("weights must be a numeric vector with one value for every
         data category, or NULL")
    }

    # which components to plot
    comp_to_plot <- colnames(ccc_mean_df)
    cat_cols <- rep(col, ncat)

    # main titles
    plot_title <- paste("Signature", comp_to_plot)

    names(plot_title) <- comp_to_plot

    for (ii in seq_along(comp_to_plot)){

      cname <- comp_to_plot[ii]

      # mean categorical distribution (sig), and credibility interval
      sig <- ccc_mean_df[,ii]
      ci <- ccc_credint[[ii]]

      # adjust categories by weights if specified (lose cred intervals though)
      if(!is.null(weights)){
        sig <- sig %*% diag(weights)
        denom <- sum(sig)
        sig <- sig/denom
        ci <- NULL # not sure how to get cred int if adjusting with weights
      }

      sig <- as.vector(sig)

      # set categories whose credibility intervals hit zero to a different colour
      cat_cols_copy <- cat_cols

      # max plotting height
      ci[is.na(ci)] <- 0
      sig[is.na(sig)] <- 0
      if(max(ci)==0){
        plottop <- max(sig)+0.05
      }else{
        plottop <- ceiling(max(ci)/0.1)*0.1+0.01
      }

      # main barplot
      b <- barplot(sig, col=cat_cols_copy, xaxt="n", ylim=c(0,plottop*1.1),
                   border=NA, names.arg=rep("", ncat), xpd=F, las=1,
                   main=plot_title[ii])

      # add credibility intervals
      if (cred_int & !is.null(ci)){
        segments(x0=b, y0=ci[1,], y1=ci[2,], col="grey30")
      }

      # add category names
      if (!is.null(cat_names)){
        mtext(cat_names, side=1, las=2, at=b, cex=cex.cat,
              family="mono", col=cat_cols)
      }
    }
  }



#' Plot the distribution of raw clusters highly similar as the component in posterior chains
#' @param components  A matrix that containing components with each row corresponding a category and each column
#'                    corresponding a component
#'
#' @param retval An object return from \code{\link{extract_ccc_cdc_from_hdp}}
#'
#' @export
plot_component_posterior_samples <- function(components,
                                             retval){
  for(i in 1:ncol(components)){
    chain <- exposures <- sequence <- NULL
    summary.cluster <- retval[[i]][["summary.chain.info"]]
    colnames(summary.cluster) <- c("chain","sample","sequence","exposures")
    cluster.name <- colnames(components)[i]
    plot.1 <- ggplot2::ggplot(data=summary.cluster, ggplot2::aes(x=sequence, y=sample, group=chain,color=chain)) +
      ggplot2::geom_point()+ggplot2::ggtitle(paste0(cluster.name," in Gibbs sample")) + ggplot2::xlab("Posterior.Sample") +  ggplot2::ylab("Chain")
    plot(plot.1)

    plot.2 <- ggplot2::ggplot(data=summary.cluster, ggplot2::aes(x=sequence, y=exposures, group=chain,color=chain)) +
      ggplot2::geom_point()+ggplot2::ggtitle(paste0("exposures of ",cluster.name," in Gibbs sample"))+ ggplot2::xlab("Posterior.Sample") +  ggplot2::ylab("Exposure")
    plot(plot.2)
  }
}

