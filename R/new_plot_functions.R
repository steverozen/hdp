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
    ccc_unlist <- do.call(cbind,ccc_0[[chain]])
    cdc_unlist <- do.call(cbind,cdc_0[[chain]])

    temp.chain <- data.frame(matrix(nrow=ncol(ccc_unlist),ncol=3))
    temp.chain[,1] <- chain
    temp.chain[,2] <- 0
    temp.chain[,3] <- c(1:nrow(temp.chain))

    cos.sims <- apply(ccc_unlist,2,function(x){lsa::cosine(x,spectrum)})
    spectrum.ccc <- cbind(spectrum.ccc,ccc_unlist[,cos.sims>cos.merge])
    temp.chain[cos.sims>cos.merge,2] <- 1
    temp.chain[cos.sims>cos.merge,4] <- colSums(ccc_unlist[,cos.sims>cos.merge])
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

#' Plot hdp signature exposure in each sample. This function returns the plot of
#' top 5 samples with the highest exposure to a component
#' @param hdpsample  A hdpSampleChain or hdpSampleMulti object including output
#'  from \code{\link{hdp_extract_components}}
#' @param retval an object return from \code{\link{extract_ccc_cdc_from_hdp}}
#' @param input.catalog input catalog for samples
#' @param col_comp Colours of each component, from 0 to the max number. If NULL,                          default colors will be used
#' @param incl_numdata_plot Logical - should an upper barplot indicating the number of
#'  data items per DP be included? (Default TRUE)
#' @param ylab_numdata Vertical axis label for numdata plot
#' @param ylab_exp Vertical exis label for exposure plot
#' @param leg.title Legend title
#' @param cex.names Expansion factor for bar labels (dpnames) in exposure plot
#' @param cex.axis Expansion factor for vertical-axis annotation
#' @param mar See ?par
#' @param oma See ?par
#' @param ... Other arguments to plot
#' @importFrom beeswarm beeswarm
#' @export
plot_sig_exposure_for_dp <- function(retval, hdpsample, input.catalog,
                                     col_comp = NULL,
                                     incl_numdata_plot=TRUE,
                                     ylab_numdata="Number of data items",
                                     ylab_exp="Component exposure",
                                     leg.title="Component", cex.names=0.6,
                                     cex.axis=0.7, mar=c(1, 4, 2, 0.5),
                                     oma=c(1.5, 1.5, 1, 1), ...){

  # input checks

  exposures <- retval$exposure

  signature <- retval$signature

  ndp <- hdpsample@chains[[1]]@hdp@numdp
  ncomp <- ncol(signature)

  dpnames <- colnames(input.catalog)

  dpindices <- (ndp-length(dpnames)+1):ndp

  if(is.null(col_comp)){
    col_comp <- grDevices::rainbow(ncol(signature), alpha = 1)
  }


  if (!is.numeric(dpindices) | any(dpindices %% 1 != 0) |
      any(dpindices < 1) | any(dpindices > ndp)) {
    stop(paste("dpindices must be integers between 1 and", ndp))
  }


  if(class(incl_numdata_plot) != "logical") {
    stop("incl_numdata_plot must be TRUE or FALSE")
  }

  # save pre-existing par conditions, and reset on exit
  par_old <- par(no.readonly=TRUE)
  on.exit(par(par_old), add=TRUE)

  # Number of data items per DP
  if (class(hdpsample) == "hdpSampleChain") {
    dps <- dp(final_hdpState(hdpsample))[dpindices]
    pps <- ppindex(final_hdpState(hdpsample))[dpindices]
  } else if (class(hdpsample) == "hdpSampleMulti") {
    dps <- dp(final_hdpState(chains(hdpsample)[[1]]))[dpindices]
    pps <- ppindex(final_hdpState(chains(hdpsample)[[1]]))[dpindices]
  }

  numdata <- sapply(dps, function(x) x@numdata)
  dp_order <- order(numdata, decreasing=TRUE)

  # if incl_numdata_plot TRUE, throw error if one DP has no data associated
  if (incl_numdata_plot & any(numdata == 0)) {
    stop("Can't have incl_numdata_plot TRUE if
         one or more dpindices have no associated data item/s")
  }

  # if different parent indices, warning
  if (length(unique(pps)) > 1) {
    warning("some dpindices have different parent nodes,
            separate plots may be better")
  }

  # which components to include in this plot
  inc <- which(rowSums(exposures, na.rm=T)>0)

  num_leg_col <- floor(sqrt(length(inc)))

  if (incl_numdata_plot){
    par(mfrow=c(2, 1), mar=mar, oma=oma, cex.axis=cex.axis, las=2)

    barplot(numdata[dp_order], main="Mutations of each tumor", col="gray", space=0, border=NA,
            names.arg="", ylab=ylab_numdata,
            legend.text=names(inc),
            args.legend=list(fill=col_comp[inc], bty="n", title=leg.title,
                             ncol=num_leg_col), ...)

    barplot(as.matrix(exposures[inc, dp_order, drop=FALSE]), space=0, col=col_comp[inc], border=NA,
            ylim=c(0, 1), names.arg=dpnames[dp_order], ylab=ylab_exp,
            cex.names=cex.names, ...)
  } else {

    par(cex.axis=cex.axis, las=2)
    # don't understand why legend.text needs rev() here and not in above case,
    # but seems to work?
    barplot(as.matrix(exposures[inc, dp_order, drop=FALSE]), space=0, col=col_comp[inc],
            border=NA, ylim=c(0, 1),
            xlim=c(0, length(dpindices) + num_leg_col + 1),
            names.arg=dpnames[dp_order],
            ylab=ylab_exp, cex.names=cex.names,
            legend.text=rev(names(inc)),
            args.legend=list(fill=col_comp[inc], bty="n", title=leg.title,
                             ncol=num_leg_col), ...)
  }

  data.exposures <- t(numdata*t(exposures))
  colnames(data.exposures) <- colnames(input.catalog)
  row.names(data.exposures) <- colnames(retval$signature)
  numdata.cutoff <- 0.5*nrow(input.catalog)
  data.exposures <- data.exposures[,which(colSums(input.catalog)>numdata.cutoff)]##exclude extremely low samples
  exposures <- exposures[,which(colSums(input.catalog)>numdata.cutoff)]
  input.catalog <- input.catalog[,which(colSums(input.catalog)>numdata.cutoff)]
  x <- barplot(rowSums(data.exposures), las=2,cex.names = 0.8) # Do not plot any axes


  Signature <- Sample <- Exposure <- Tumor <- NULL

  df <- reshape2::melt(data.exposures)
  colnames(df) <- c("Signature","Sample","Exposure")

  if(any(grepl("::",colnames(data.exposures)))){

    df$Tumor <- apply(df,1,function(x){
      x["Tumor"] <- unlist(strsplit(x["Sample"],"::"))[1]
    })
  }else{
    df$Tumor <- "No specific tumor type"
  }

  ##for each signature, plotting the top 5 tumors
  for(i in 1:nrow(data.exposures)){
    dp_order_sig <- order(exposures[i,], decreasing=TRUE)
    this.par <- par(mfrow=c(2, 1), mar=mar, oma=oma, cex.axis=cex.axis, las=1)

    this.col <- rep("grey",length(col_comp))
    this.col[i] <- col_comp[i]


    sig <- row.names(data.exposures)[i]

    plot(exposures[i,]~data.exposures[i,],xlab="Exposure",ylab="Prop.Exposure",main=sig,pch=16)

    #sig.df <- df[df$Signature==sig,]

    on.exit(par(this.par))

    #sig.df$exp.prop <- as.numeric(exposures[i,])
    #beeswarm::beeswarm(Exposure~Tumor,data=sig.df,method="swarm",col=1:3, pch=19, cex=.75)
    #ggplot2::ggplot(data=sig.df, aes(x=Tumor, y=exp.prop)) +
    #  ggplot2::geom_bar(stat="identity")+ ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))

    old.par <- par(mfrow = c(6, 1), mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
    on.exit(par(old.par))

    ICAMS::PlotCatalog(ICAMS::as.catalog(signature[,i,drop=FALSE],infer.rownames = T,catalog.type = "counts.signature"))

    this.catalog <- input.catalog[,dp_order_sig[1:5], drop=FALSE]
    this.prop <- exposures[i,dp_order_sig[1:5]]
    this.prop <- round(this.prop,3)
    colnames(this.catalog) <- paste0(colnames(this.catalog),"(",this.prop,")")
    for (j in 1:5) {

      ICAMS::PlotCatalog(this.catalog[,j, drop=FALSE])
    }


  }

}

#' Plot the distribution of raw clusters highly similar as the component in posterior chains
#' @param components  A matrix that containing components with each row corresponding a category and each column
#'                    corresponding a component
#'
#' @param retval An object return from \code{\link{extract_ccc_cdc_from_hdp}}
#'
plot_component_posterior_samples <- function(components,
                                             retval){
  for(i in 1:ncol(signatures)){
    chain <- signatures <- exposures <- NULL
    summary.cluster <- retval[[i]][["summary.chain.info"]]
    colnames(summary.cluster) <- c("chain","sample","sequence","exposures")
    cluster.name <- colnames(signatures)[i]
    plot.1 <- ggplot2::ggplot(data=summary.cluster, ggplot2::aes(x=sample, y=sequence, group=chain,color=chain)) +
      ggplot2::geom_point()+ggplot2::ggtitle(paste0(cluster.name," in Gibbs sample")) + ggplot2::xlab("Posterior.Sample") +  ggplot2::ylab("Chain")
    plot(plot.1)

    plot.2 <- ggplot2::ggplot(data=summary.cluster, ggplot2::aes(x=sample, y=exposures, group=chain,color=chain)) +
      ggplot2::geom_point()+ggplot2::ggtitle(paste0("exposures of ",cluster.name," in Gibbs sample"))+ ggplot2::xlab("Posterior.Sample") +  ggplot2::ylab("Exposure")
    plot(plot.2)
  }
}
