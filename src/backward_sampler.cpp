#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// use the cholesky decomposition to draw a value from a MV Gaussian distribution
//' @export
// [[Rcpp::export]]
arma::colvec mvrnorm(const arma::colvec mean, const arma::mat Sigma) {
  const int n = mean.n_elem;
  arma::colvec z(n);
  for (int i = 0; i < n; ++i) {
    z.at(i) = R::rnorm(0.0, 1.0);
  }
  
  // Calculate the cholesky decomposition of Sigma
  arma::mat root = arma::chol(Sigma);
  
  // draw from a multivariate normal
  arma::colvec x = root * z + mean;
  return x;
}

//' @export
// [[Rcpp::export]]
arma::colvec mvrnormsvd(
  const arma::colvec mean,
  const arma::mat Sigma) {
  const int n = mean.n_elem;
  arma::colvec z(n);
  for (int i = 0; i < n; ++i) {
    z.at(i) = R::rnorm(0.0, 1.0);
  }
  
  arma::mat eigenvectors;
  arma::colvec eigenvalues;
  arma::eig_sym(eigenvalues, eigenvectors, Sigma);
  
  return mean + (eigenvectors * arma::diagmat(arma::sqrt(eigenvalues)) * z);
}

arma::colvec samplingStep(
    const arma::mat g,
    const arma::colvec mt,
    const arma::mat ct,
    const arma::colvec at1,
    const arma::mat rt1,
    const arma::colvec theta,
    const arma::mat w) {
  
  // more efficient than inverting rt, equivalent to C * G.t * inv(R)
  const arma::mat cgrinv = solve(rt1.t(), g * ct.t()).t();
  
  // calculate the updated mean
  const arma::mat mean = mt + cgrinv * (theta - at1);
  
  // calculate the updated covariance using Joseph Form
  int n = mt.n_rows;
  const arma::mat identity = arma::eye(n, n);
  const arma::mat diff = identity - cgrinv * g;
  const arma::mat covariance = diff * ct * diff.t() + cgrinv * w * cgrinv.t();
  
  // ensure symmetry
  const arma::mat r = (covariance + covariance.t()) / 2.0;
  
  return mvrnormsvd(mean, r);
}

//' @export
// [[Rcpp::export]]
arma::mat dlm_backward_sample(
    const arma::mat g,
    const arma::mat w,
    const arma::mat mts,
    const arma::cube cts,
    const arma::mat ats,
    const arma::cube rts) {
  
  arma::colvec mt0 = mts.col(0);
  int n = mts.n_cols - 1;
  int p = mt0.n_rows;
  
  arma::mat theta(p, n + 1, arma::fill::none);
  theta.col(n) = mvrnorm(mts.col(n), cts.slice(n));

  for (int t = n - 1; t >= 0; --t) {
    arma::colvec mt = mts.col(t);
    arma::mat ct = cts.slice(t);
    arma::colvec at = ats.col(t);
    arma::mat rt = rts.slice(t);

    theta.col(t) = samplingStep(g, mt, ct, at, rt, theta.col(t + 1), w);
  }

  return theta;
}
