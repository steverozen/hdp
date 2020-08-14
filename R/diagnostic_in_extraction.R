#' Diagnostic plotting inside of hdp_merge_and_extract_components function.
#' This function generates details of the raw clusters in hdp.0
#'
#' @param clust_hdp0_ccc An object in \code{\link{hdp_merge_and_extract_components}}.
#' @param ccc clust_categ_counts from \code{\link{hdp_merge_and_extract_components}}.
#' @param cdc clust_dp_counts from \code{\link{hdp_merge_and_extract_components}}.
#' @param diagnostic.folder A directory where details for hdp.0 are plotted.
#' @param ncat Number of categories. An internal object from \code{\link{hdp_merge_and_extract_components}}.
#' @param nsamp Number of posterior samples. An internal object from \code{\link{hdp_merge_and_extract_components}}.
#' @param nch Number of posterior chains. An internal object from \code{\link{hdp_merge_and_extract_components}}.
#'
#'
#'
#' @return The plots of presence of a raw cluster in each chain.
#' @seealso \code{\link{hdp_merge_and_extract_components}}
#'
#' @export


diagnostic_in_extraction <- function(clust_hdp0_ccc,
                                     ncat,
                                     nsamp,
                                     nch,
                                     ccc,
                                     cdc,
                                     diagnostic.folder){

  chain <- exposures <- nsampchain <- NULL

  nsampchain <- nsamp/nch
  for(i in 1:ncol(clust_hdp0_ccc)){
    summary.matrix <- data.frame(colnames(clust_hdp0_ccc))
    cluster.pattern <- clust_hdp0_ccc[,i]
    cluster.name <- colnames(clust_hdp0_ccc)[i]


    ccc_cluster <- which(colnames(ccc[[1]]) == unlist(strsplit(cluster.name,"[_]",perl=T))[3])


    summary.cluster <- data.frame(matrix(ncol=0,nrow=nsamp))

    individual.catalog <- data.frame(matrix(nrow=ncat,ncol=0))


    for(index.i in 1:nch){

      for(index.j in 1:nsampchain){

        index <- (index.i-1)*nsampchain+index.j

        temp.matrix <- ccc[[index]]
        temp.cdc.matrix <- cdc[[index]]

        summary.cluster$chain[index] <- paste("chain.",index.i,sep="")

        summary.cluster$sequence[index] <- 0
        summary.cluster$exposures[index] <- 0
        summary.cluster$sample[index] <- index.j


        if(sum(temp.matrix[,ccc_cluster])>0){
          cosine <- lsa::cosine(temp.matrix[,ccc_cluster],cluster.pattern)
          individual.catalog <- cbind(individual.catalog,temp.matrix[,ccc_cluster])
          if(cosine>0.95){
            summary.cluster$sequence[index] <- index.i
            summary.cluster$exposures[index] <- sum(temp.cdc.matrix[,ccc_cluster])

          }
        }


      }

    }

    dir.create(file.path(diagnostic.folder,cluster.name), recursive = T)

    row.names(individual.catalog) <- ICAMS::catalog.row.order$SBS96
    if(ncol(individual.catalog)>0){
      individual.catalog.catalog <- ICAMS::as.catalog(individual.catalog,catalog.type = "counts")
      ICAMS::PlotCatalogToPdf(individual.catalog.catalog,
                              file.path(diagnostic.folder,cluster.name,"/individual.raw.cluster.pdf"))

    }
    grDevices::pdf(file.path(diagnostic.folder,cluster.name,"/cluster.in.Gibbs.sample.pdf"))
    plot.1 <- ggplot2::ggplot(data=summary.cluster, ggplot2::aes(x=sample, y=sequence, group=chain,color=chain)) +
      ggplot2::geom_point()+ggplot2::ggtitle(paste0(cluster.name," in Gibbs sample")) + ggplot2::xlab("Posterior.Sample") +  ggplot2::ylab("Chain")
    plot(plot.1)
    grDevices::dev.off()

    grDevices::pdf(file.path(diagnostic.folder,cluster.name,"/cluster.exposure.in.Gibbs.sample.pdf"))
    plot.2 <- ggplot2::ggplot(data=summary.cluster, ggplot2::aes(x=sample, y=exposures, group=chain,color=chain)) +
      ggplot2::geom_point()+ggplot2::ggtitle(paste0("exposures of ",cluster.name," in Gibbs sample"))+ ggplot2::xlab("Posterior.Sample") +  ggplot2::ylab("Exposure")
    plot(plot.2)
    grDevices::dev.off()

  }

}
