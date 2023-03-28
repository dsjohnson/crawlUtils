#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector kde_estimate(NumericMatrix grid,
                           NumericMatrix points,
                           NumericMatrix B,
                           bool norm = true) {

  int N = grid.nrow();
  int n = points.nrow();
  double kern_sum;
  NumericVector out(N);
  NumericVector tmp(N);
  NumericVector y;
  NumericVector x;
  NumericVector d;
  for (int j = 0; j < n; j++) {
    NumericVector x = points.row(j);
    for (int i = 0; i < N; i++) {
      y = grid.row(i);
      d = y-x;
      tmp(i) = exp(-0.5*( d(0)*d(0)*B(0,0) + 2*d(0)*d(1)*B(0,1) + d(1)*d(1)*B(1,1) ));
    }
    if(norm){
      kern_sum = sum(tmp);
      if(kern_sum>0) tmp = tmp/kern_sum;
    }
    out += tmp;
    checkUserInterrupt();
  }
  out = out/sum(out);
  return out;
}
