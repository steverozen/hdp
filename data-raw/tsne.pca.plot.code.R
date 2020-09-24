#' A pca plot for signatures and tumors
#'
#' @param exposure An exposure matrix
#' @param colorcategory A list of category of samples for labelling tsne plot.
#'                Length is equal as number of samples
#' @importFrom graphics abline axis barplot matplot mtext par plot points segments text
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes element_text geom_point theme
# @examples
#' @export
#' @rdname plotcomp
plot_pca_sigs_tumortype <- function(exposure,
                                    colorcategory = NULL){
  ColorCat <- x <- y <- NULL
  pca.out <- stats::prcomp(t(exposure), scale. = TRUE)
  pca.plot <- data.frame(x=pca.out$x[,1],y=pca.out$x[,2],
                         tumorname = colnames(exposure))

  if(!is.null(colorcategory)){
    pca.plot$ColorCat <- colorcategory
  }else{
    pca.plot$ColorCat <- apply(pca.plot,1,function(x){
      x["ColorCat"] <- unlist(strsplit(x["tumorname"],"::"))[1]
    })
  }

  plot <- ggplot2::ggplot(pca.plot) + ggplot2::geom_point(ggplot2::aes(x=x, y=y, color=ColorCat))+
    ggplot2::theme(legend.text= ggplot2::element_text(size=8))+ggplot2::ggtitle("All Signatures")

  plot(plot)

  for(sig in row.names(exposure)){
    pca.plot$sig <- exposure[sig,]
    plot <- ggplot2::ggplot(pca.plot) + ggplot2::geom_point(ggplot2::aes(x=x, y=y, color=sig))+
      ggplot2::theme(legend.text= ggplot2::element_text(size=8))+ggplot2::ggtitle(sig)+
      ggplot2::scale_color_gradient(low = "yellow", high = "red")
    plot(plot)

  }

}



#' A tsne plot for signatures and tumors
#'
#' @param exposure An exposure matrix
#' @param colorcategory A list of category of samples for labelling tsne plot.
#'                Length is equal as number of samples
#' @importFrom graphics abline axis barplot matplot mtext par plot points segments text
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes element_text geom_point theme scale_color_gradient
#' @importFrom Rtsne Rtsne
# @examples
#' @export
#' @rdname plotcomp
plot_tsne_sigs_tumortype <- function(exposure,
                                     colorcategory = NULL){
  ColorCat <- x <- y <- NULL
  set.seed(44) # for reproducibility
  tsne.out<- Rtsne::Rtsne(t(exposure),check_duplicates = FALSE)
  tsne.plot <- data.frame(x = tsne.out$Y[,1], y = tsne.out$Y[,2], tumorname = colnames(exposure))

  if(!is.null(colorcategory)){
    tsne.plot$ColorCat <- colorcategory
  }else{
    tsne.plot$ColorCat <- apply(tsne.plot,1,function(x){
      x["ColorCat"] <- unlist(strsplit(x["tumorname"],"::"))[1]
    })
  }


  plot <- ggplot2::ggplot(tsne.plot) + ggplot2::geom_point(ggplot2::aes(x=x, y=y, color=ColorCat))+
    ggplot2::theme(legend.text= ggplot2::element_text(size=8))+ggplot2::ggtitle("All Signatures")

  plot(plot)

  for(sig in row.names(exposure)){
    tsne.plot$sig <- exposure[sig,]
    plot <- ggplot2::ggplot(tsne.plot) + ggplot2::geom_point(ggplot2::aes(x=x, y=y, color=sig))+
      ggplot2::theme(legend.text= ggplot2::element_text(size=8))+ggplot2::ggtitle(sig)+
      ggplot2::scale_color_gradient(low = "yellow", high = "red")
    plot(plot)

  }

  ##ToDO a pca for every signature exposure


}

