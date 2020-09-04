// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix cosCpp(NumericMatrix Xr) {
  int n = Xr.nrow(), k = Xr.ncol();
  mat X(Xr.begin(), n, k, false); // reuses memory and avoids extra copy
  mat Y = trans(X) * X; // matrix product
  mat res = Y / (sqrt(diagvec(Y)) * trans(sqrt(diagvec(Y))));
  return wrap(res);
}
