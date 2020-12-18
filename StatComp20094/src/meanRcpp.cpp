#include <Rcpp.h>
using namespace Rcpp;

//' @title Mean value
//' @description Calculate the mean value using Rcpp
//' @param x a numericVector of samples
//' @return mean value \code{m}
//' @examples
//' \dontrun{
//' x<-1:100
//' meanRcpp(x)
//' }
//' @export
// [[Rcpp::export]]
double meanRcpp(NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}

