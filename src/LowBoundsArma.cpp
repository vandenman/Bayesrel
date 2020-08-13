//' @importFrom Rcpp evalCpp
//' @useDynLib Bayesrel, .registration=TRUE
//'
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//[[Rcpp::export]]
double alphaArma(arma::mat X) {
    double k = X.n_cols;
	double tr = arma::trace(X);
	double out = k/(k-1) * (1 - (tr/(arma::accu(X))));
	return out;
}

//[[Rcpp::export]]
double l2Arma(arma::mat X) {
  double k = X.n_cols;
  mat X0 = X;
  X0.diag().zeros();
  double out = (accu(X0) + sqrt(k/(k-1) * accu(square(X0)))) / accu(X);
  return out;
}

//[[Rcpp::export]]
double l6Arma(arma::mat X) {
  // correlation matrix from covariance matrix:
  vec sds = 1/sqrt(X.diag());
  mat Xcor = diagmat(sds) * X * diagmat(sds);
  Xcor.diag().ones();
  mat XCorInv = inv_sympd(Xcor);
  vec smc = 1 - 1 / XCorInv.diag();
  double out = 1 - accu(1 - smc) / accu(Xcor);
  return out;
}
