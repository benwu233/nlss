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
NumericMatrix Simu_DLSS(NumericMatrix S,NumericMatrix A, int K){
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
          if( S(h,j)==(k+1) ){
            prob[k] += A(i,h);
          }
        }
      }
      out(i,j) = as<double>(Rcpp::sample(K,1,false,as<sugar::probs_t>(prob)));
    }
  }

  return out;
}


// [[Rcpp::export]]
NumericMatrix Simu_data_bi(NumericMatrix X, NumericMatrix mu0, NumericMatrix sd0, NumericVector th){
  Rcpp::Dimension X_dim = X.attr("dim");
  int n = X_dim[0];
  int p = X_dim[1];
  IntegerVector n1(1,0);
  IntegerVector n2(1,0);
  int tag1 = 0;
  int tag2 = 0;

  NumericMatrix out(n,p);

  NumericVector mu(1,0.0);
  NumericVector sd(1,0.0);
  NumericVector a(1,0.0);
  NumericVector b(1,0.0);
  NumericVector mu1(1,0.0);
  NumericVector sd1(1,0.0);
  NumericVector a1(1,0.0);
  NumericVector b1(1,0.0);

  Rcpp::Environment TGaussian_env("package:truncnorm");
  Rcpp::Function rtruncnorm_r = TGaussian_env["rtruncnorm"];
  for(int i = 0; i < n; i++){
    mu[0] = mu0(i,0);
    sd[0] = sd0(i,0);
    a[0] = -100;
    b[0] = th[i];

    mu1[0] = mu0(i,1);
    sd1[0] = sd0(i,1);
    a1[0] = th[i];
    b1[0] = 100;

    n1[0] = 0;
    n2[0] = 0;
    for(int j = 0; j < p; j++){
      if(X(i,j)==1){
        n1[0] ++;
      }
      else if(X(i,j)==2){
        n2[0] ++ ;
      }
    }

    NumericVector candidate1(n1[0],0.0);
    NumericVector candidate2(n2[0],0.0);

    candidate1 = rtruncnorm_r(n1,a,b,mu,sd);
    candidate2 = rtruncnorm_r(n2,a1,b1,mu1,sd1);

    tag1 = 0;
    tag2 = 0;
    for(int j = 0; j < p; j++){
      if(X(i,j)==1){
        out(i,j) = candidate1[tag1];
        tag1 ++ ;
      }
      else if(X(i,j)==2){
        out(i,j) = candidate2[tag2];
        tag2 ++;
      }
    }
  }


  return out;

}


// [[Rcpp::export]]
double DLSS_Deviance_c(NumericMatrix X,NumericMatrix S,NumericMatrix A, int K){
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];

  double Dev = 0;
  double prob = 0;

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
          Dev += log(prob + 1e-20);
        }
      }
    }
  }

  Dev = Dev * (-2);

  return Dev;

}



// [[Rcpp::export]]
double DLSS_logLik_noisefree(NumericMatrix X,NumericMatrix S,NumericMatrix A, int K){
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];

  double out = 0;
  double prob = 0;

  for(int k=0; k<K; k++){
    for(int j=0; j<p; j++){
      for(int i=0; i<n; i++){
        if( X(i,j)==(k+1) ){
          prob = 0;
          for(int h=0 ; h<q; h++){
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

// [[Rcpp::export]]
double NLSS_logLik_noise(NumericMatrix X,NumericMatrix A,NumericMatrix beta,IntegerVector group, int K, int G){
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];

  double out = 0;
  double prob = 0;
  int k = 0;
  int g = 0;

  for(int j=0; j<p; j++){
    for(int i=0; i<n; i++){

      prob = A(i,q-1)/K;

      for(int h=0 ; h<(q-1); h++){
        g = group[j];
        k = (int)(X(i,j) - 1 + 1e-20);
        prob += A(i,h) * beta(h+g*(q-1),k);
      }
      out += log(prob);
    }
  }

  return out;
}



// [[Rcpp::export]]
double NLSS_logLik_noise_group(NumericMatrix X,NumericMatrix A,NumericMatrix beta,IntegerVector group, int g0, int K, int G){
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];

  double out = 0;
  double prob = 0;
  int k = 0;
  int g = 0;

  for(int j=0; j<p; j++){

    g = group[j];
    if(g==g0){
    for(int i=0; i<n; i++){

      prob = A(i,q-1)/K;

      for(int h=0 ; h<(q-1); h++){
          k = (int)(X(i,j) - 1 + 1e-20);
          prob += A(i,h) * beta(h+g*(q-1),k);
      }
      if(prob > 0){
        out += log(prob);
      }

    }
    }

  }

  return out;
}

// [[Rcpp::export]]
double DLSS_logLik_noise_1(NumericMatrix X,NumericMatrix S,NumericMatrix A,NumericMatrix beta,IntegerVector group, int K, int G){
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A.attr("dim");
  Rcpp::Dimension S_dim = S.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = S_dim[0];

  double out = 0;
  double prob = 0;
  double count = 0;
  int k = 0;
  int g = 0;

  for(int j=0; j<p; j++){
    for(int i=0; i<n; i++){

      k = (int)(X(i,j) - 1 + 1e-20);
      prob = A(i,q+k)/K;

      for(int h=0 ; h<q; h++){

        g = group[j];
        prob += A(i,h) * beta(h+g*q,k);

      }


      out += log(prob + 1e-20);
    }

  }

  return out;
}


// [[Rcpp::export]]
double DLSS_logLik_noise0_group(NumericMatrix X,NumericMatrix S,NumericMatrix A,IntegerVector group, int g0, int K, int G){
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];

  double out = 0;
  double prob = 0;
  double count = 0;
  int g = 0;

  for(int k=0; k<K; k++){

    for(int j=0; j<p; j++){

      g=group[j];

      if(g==g0){

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

  }

  /*
  for(int l = 0; l < (q-1); l++){
    for(int k = 0; k < K; k++){
      count  = 0;
      for(int j = 0; j < p; j++){
        if(S(l,j)==(k+1)){
          count ++;
        }
      }

      out += count * log(beta(l,k)+ 1e-20);
    }
  }

*/

  return out;

}



// [[Rcpp::export]]
double DLSS_logLik_noise0(NumericMatrix X,NumericMatrix S,NumericMatrix A,NumericMatrix beta, int K){
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

  /*
   for(int l = 0; l < (q-1); l++){
   for(int k = 0; k < K; k++){
   count  = 0;
   for(int j = 0; j < p; j++){
   if(S(l,j)==(k+1)){
   count ++;
   }
   }

   out += count * log(beta(l,k)+ 1e-20);
   }
   }

   */

  return out;

}
