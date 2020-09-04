
#' Identify individual clusters from ccc.
#'
#' @param main.index the index of cluster in a given ccc
#'
#' @param ccc a list clust_categ_counts matrices
#'
#' @param cdc a list clust_dp_counts matrices
#'
#' @param disregard index of parent nodes to be excluded
#'
#' @export

identify_individual_cluster <- function(main.index,
                                        ccc,
                                        cdc,
                                        disregard){


  cluster.name <- paste0("potential.hdp.",main.index)
  message("doing ",cluster.name)
  # dir.create(file.path(diagnostic.folder,cluster.name), recursive = T)
  compii <- sapply(ccc, function(x) x[,main.index])
  clust_label <- 1:ncol(compii)
  colnames(compii) <- clust_label

  clust_label <- generate_label_high_cossim(clust_label = clust_label,
                                            matrix = compii,
                                            cos.sim = 0.90)

  merged.compii <- data.frame(merge_cols(as.matrix(compii),clust_label))

  merged.compii <- merged.compii[,colSums(merged.compii)>0]

  merged.compii <- data.frame(merged.compii)

  if(ncol(merged.compii)>1){ ##another step of merge
    clust_label <- 1:ncol(merged.compii)
    colnames(merged.compii) <- clust_label
    clust_label <- generate_label_high_cossim(clust_label = clust_label,
                                              matrix = merged.compii,
                                              cos.sim =0.90)
    merged.compii <- data.frame(merge_cols(as.matrix(merged.compii),clust_label))
  }

  ##temporary selection. not sure if this is a good number though




  sum.list <-{}
  merged.compii <- merged.compii[,colSums(merged.compii)>(1*length(ccc))]
  merged.compii <- data.frame(merged.compii)


  if(ncol(merged.compii)>0){
    message("found ",ncol(merged.compii)," potential signatures in ",cluster.name)
    for(this.clust in 1:ncol(merged.compii)){
      colnames(merged.compii)[this.clust] <- this.clust
      target.pattern <- merged.compii[,this.clust]

      if(sum(target.pattern)>0){

        exp.df <- data.frame(matrix(nrow=nrow(cdc[[1]][-disregard,]),ncol=0))
        ccc.agg <- stats.list <- exp.list  <- {}

        for(this.samp in 1:length(ccc)){

          res <- apply(ccc[[this.samp]]+0.001, 2, lsa::cosine, y=target.pattern)

          if(any(res>0.95)){
            ccc.agg <- cbind(ccc.agg,ccc[[this.samp]][,which(res>0.95)])

            exp.df <- cbind(exp.df,cdc[[this.samp]][-disregard,which(res>0.95)])
            stats.list <- c(stats.list,this.samp)
            exp.list <- c(exp.list,sum(ccc[[this.samp]][,which(res>0.95)]))

          }

        }

      }
      if(length(stats.list)>0){


        test <- list(counts.spec = target.pattern,
                     matched.samp = length(stats.list),
                     mean.exp = mean(exp.list),
                     sd.exp = sd(exp.list),
                     exp.df = exp.df,
                     ccc.agg = ccc.agg)
      }

      sum.list <- append(sum.list, list(test))
      names(sum.list)[this.clust] <- this.clust

    }

  }else{
    message("found 0 potential signatures in ",cluster.name)
  }

  return(sum.list)

}



#'label_based_on_cosine_similarity
#'
#'
#'@param clust_label column names of the input matrix
#'@param matrix the matrix needs to be measured
#'@param cos.sim cosine similarity threshold
#'
#'@keywords  internal
generate_label_high_cossim <- function(clust_label,
                                       matrix,
                                       cos.sim){
  clust_cos <- cosCpp(as.matrix(matrix))
  clust_same <- (clust_cos > cos.sim & lower.tri(clust_cos))
  same <- which(clust_same, arr.ind=TRUE) # merge these columns
  if (length(same)>0){
    for (index in 1:nrow(same)){
      clust_label[same[index, 1]] <- clust_label[same[index, 2]]
    }
  }
  return(clust_label)
}
