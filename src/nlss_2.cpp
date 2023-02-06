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
NumericMatrix simNLSS(NumericMatrix S,NumericMatrix A, int K){
  Rcpp::Dimension S_dim = S.attr("dim");
  Rcpp::Dimension A_dim = A.attr("dim");

  int n = A_dim[0];
  int p = S_dim[1];
  int q = A_dim[1];
  NumericVector prob(K,0.0);

  NumericMatrix out(n,p);

  for(int i=0; i<n; i++){
    for(int j=0; j<p; j++){
      for(int k=0; k<K; k++){
        prob[k] = A(i,q-1)/K;
        for(int h=0 ; h<(q-1); h++){
          if( S(h,j)==(k) ){
            prob[k] += A(i,h);
          }
        }
      }
      out(i,j) = as<double>(Rcpp::sample(K,1,false,as<sugar::probs_t>(prob))) - 1;
    }
  }

  return out;
}


// [[Rcpp::export]]
double NLSS_logLik_noise0(NumericMatrix X,NumericMatrix S,NumericMatrix A, int K){
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];

  double out = 0;
  double prob = 0;
  double count = 0;

  for(int k=0; k<K; k++){

    for(int j=0; j<p; j++){

      for(int i=0; i<n; i++){
        if( X(i,j)==(k+1) ){
          prob = A(i,q-1)/K;
          for(int h=0 ; h<(q-1); h++){

            if( S(h,j)==(k+1) ){
              prob += A(i,h);
            }

          }
          out += log(prob + 1e-20);
        }
      }
    }
  }

  return out;

}
