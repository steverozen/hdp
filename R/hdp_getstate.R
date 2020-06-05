# extract key info from this hdp iteration (list class)
# numclass = number of clusters
# classqq = category vs cluster counts overall (matrix) category is the mutation type (e.g. ACT > AAT), so each column is like a proto signature profile
# classnd = dp vs cluster counts (matrix) rows are processes (samples plus other nodes), columns are clusters (proto-signatures)
# alpha = conparam values (vector)

hdp_getstate <- function(hdp){
  hdpstate <- list()
  hdpstate$numclass <- hdp$base$numclass
  hdpstate$classqq  <- hdp$base$classqq
  hdpstate$classnd  <- t(sapply(hdp$dp, function(x) x$classnd))
  hdpstate$alpha    <- sapply(hdp$conparam, function(x) x$alpha)
  return(hdpstate)
}
