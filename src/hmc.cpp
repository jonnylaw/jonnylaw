#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector hmc_leapfrog_cpp(
  Function gradient,
  NumericMatrix ys,
  NumericVector qp,
  double stepSize) {

  // unasign values
  int d = qp.size() / 2;
  NumericVector momentum(d);
  NumericVector position(d);
  for (int i = 0; i < d; i++) {
    momentum(i) = qp(d + i);
    position(i) = qp(i);
  }

  NumericVector momentum1 = momentum + as<NumericVector>(gradient(ys, position)) * 0.5 * stepSize;
  NumericVector position1 = position + stepSize * momentum1;
  NumericVector momentum2 = momentum + as<NumericVector>(gradient(ys, position1)) * 0.5 * stepSize;

  NumericVector newqp(2 * d);
  for (int i = 0; i < d; i++) {
    newqp(i) = position1(i);
    newqp(d + i) = momentum2(i);
  }
  return newqp;
}
