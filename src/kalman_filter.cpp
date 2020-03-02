#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
List advState(
    const arma::mat g, 
    const arma::colvec mt,
    const arma::mat ct,
    const arma::mat w) {
  
  const arma::colvec at = g * mt;
  const arma::mat rt = g * ct * g.t() + w;
  
  return List::create(at, rt);
}

//' @export
// [[Rcpp::export]]
List oneStep(
    const arma::mat ff,
    const arma::mat v,
    const arma::colvec at,
    const arma::mat rt) {

  const arma::mat fft = ff.t();
  const arma::colvec ft = ff * at;
  const arma::mat qt = ff * rt * fft + v;

  return List::create(ft, qt);
}

//' Update the state using Joseph Form Update given the newly observed data.
//' @export
// [[Rcpp::export]]
List updateState(
    const arma::mat ft,
    const arma::mat at,
    const arma::mat rt,
    const arma::colvec predicted,
    const arma::mat predcov,
    const arma::colvec y,
    const arma::mat v) {

  const arma::colvec residual = y - predicted;

  const arma::mat kalmanGain = solve(predcov.t(), (ft * rt.t())).t();
  const arma::colvec mt1 = at + kalmanGain * residual;
  int n = mt1.n_rows;
    
  const arma::mat identity = arma::eye(n, n);

  const arma::mat diff = (identity - kalmanGain * ft);
  const arma::mat covariance = diff * rt * diff.t() + kalmanGain * v * kalmanGain.t();

  return List::create(mt1, covariance);
}

//' @export
// [[Rcpp::export]]
arma::mat missing_f(const arma::mat f, const arma::colvec y) {
  return f.rows(arma::find_finite(y));
}

//' @export
// [[Rcpp::export]]
arma::mat missing_v(const arma::mat v, const arma::colvec y) {
  return v.submat(arma::find_finite(y), arma::find_finite(y));
}

//' @export
// [[Rcpp::export]]
arma::colvec flatten(arma::colvec X) {
  return X.elem(arma::find_finite(X));
}

//' @export
// [[Rcpp::export]]
List dlm_kalman_step(
    const arma::colvec yo,
    const arma::mat f,
    const arma::mat g,
    const arma::mat v,
    const arma::mat w,
    const arma::colvec mt,
    const arma::mat ct) {
  arma::colvec y = flatten(yo);
  List atrt = advState(g, mt, ct, w); //advance the state
  arma::colvec at = atrt[0];
  arma::mat rt = atrt[1];
  
  List ftqt = oneStep(f, v, at, rt); // calculate the one-step prediction
  arma::colvec ft = ftqt[0];
  arma::mat qt = ftqt[1];

  if (yo.is_finite()) { // if there is no non-missing data
    List mtct = updateState(f, at, rt, ft, qt, yo, v);
    arma::colvec mt1 = mtct[0];
    arma::mat ct1 = mtct[1];
    
    return List::create(mt1, ct1, at, rt);
  } else if (y.n_rows == 0) { // if the data is totally missing
    return List::create(at, rt, at, rt);
  } else { // if there is partially missing data
    arma::mat vm = missing_v(v, y);
    arma::mat fm = missing_f(f, y);
    List predlist = oneStep(fm, vm, at, rt);
    arma::colvec predicted = predlist[0];
    arma::mat predcov = predlist[1];

    List mtct = updateState(fm, at, rt, predicted, predcov, y, vm);
    arma::colvec mt1 = mtct[0];
    arma::mat ct1 = mtct[1];
    return List::create(mt1, ct1, at, rt);
  }
}

//' @export
// [[Rcpp::export]]
List dlm_kalman_filter(
    const arma::mat ys,
    const arma::mat f,
    const arma::mat g,
    const arma::mat v,
    const arma::mat w,
    const arma::colvec m0,
    const arma::mat c0) {
  int n = ys.n_cols;
  int d = m0.n_rows;
  arma::mat mt(d, n + 1);
  arma::cube ct(d, d, n + 1);
  arma::mat at(d, n);
  arma::cube rt(d, d, n);

  // initialise the state
  mt.col(0) = m0;
  ct.slice(0) = c0;

  for (int t = 0; t < n; ++t) {
    List mtct = dlm_kalman_step(ys.col(t), f, g, v, w, mt.col(t), ct.slice(t));
    arma::colvec mt1 = mtct[0];
    arma::mat ct1 = mtct[1];
    arma::colvec at1 = mtct[2];
    arma::mat rt1 = mtct[3];
    mt.col(t + 1) = mt1;
    ct.slice(t + 1) = ct1;
    at.col(t) = at1;
    rt.slice(t) = rt1;
  }
  
  List L = List::create(Rcpp::Named("mt")=mt, 
                        Rcpp::Named("ct")=ct, 
                        Rcpp::Named("at")=at, 
                        Rcpp::Named("rt")=rt);
  
  L.attr("class") = "filtered";
  
  return L;
}