#include <Rcpp.h>
#include <RcppParallel.h>
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


// [[Rcpp::depends(RcppParallel)]]
using namespace RcppParallel;

struct UpdateY_n : public Worker {
  RVector<double> Y;

  RVector<double> X;
  RVector<double> A;
  RVector<double> S;
  RVector<double> seed;

  int q;
  int p;
  int n;
  int K;

  UpdateY_n(NumericVector Y, NumericVector X, NumericVector A, NumericVector S,
          NumericVector seed, int q, int p, int n, int K)
    : Y(Y), X(X), A(A), S(S), seed(seed),
      q(q), p(p), n(n), K(K) {}


  void operator() (std::size_t begin, std::size_t end) {

    double prob_Y[q];
    double sum_prob_Y = 0;
    int ind = 0;

    double tmp = 0;
    double acm = 0;

    for(std::size_t j=begin;j<end;j++){
      for(std::size_t i=0;i<n; i++){

        sum_prob_Y = 0;
        ind = i+n*j;
        for(std::size_t l = 0; l<(q-1); l++){
          if(X[ind]==S[l+(q-1)*j]){
            prob_Y[l] = A[i+n*l];
          }
          else{
            prob_Y[l] = 0;
          }
          sum_prob_Y += prob_Y[l];
        }
        prob_Y[q-1] = A[i+n*(q-1)];
        sum_prob_Y += prob_Y[q-1]/K;

        tmp = seed[ind] * sum_prob_Y;
        acm = 0;

        for(std::size_t l = 0; l<q; l++){
          acm += prob_Y[l];
          if(acm > tmp){
            Y[i+n*j] = l + 1;
            break;
          }
        }
      }
    }
  }
};

inline double DLSS_compute_loglik_pll_n(RVector<double> X,RVector<double> A, RVector<double> beta,
                                      int j, int K, int n, int q, double *S_j){

  j = j - 1;
  double loglik = 0.0;
  double temp = 0.0;
  int S_jm = 0;

  for(std::size_t kk=0;kk<K; kk++){
    for(std::size_t i=0;i<n; i++){
      if(X[i+n*j]==kk+1){
        temp = 0.0;
        for(std::size_t m=0;m<(q-1);m++){
          if(S_j[m]==kk+1){
            temp += A[i+n*m];
          }
        }
        temp += A[i+n*(q-1)]/K;
        loglik += log(temp+1e-20);
      }
    }
  }
  for(std::size_t m=0;m<(q-1); m++){
    S_jm = (int)(S_j[m]+1e-20);
    loglik += log(beta[m*K + S_jm - 1]);
  }
  return loglik;
}

struct UpdateS_n : public Worker {
  RVector<double> S;
  RVector<double> loglik;

  RVector<double> X;
  RVector<double> A;

  RVector<double> beta;

  RVector<double> u;
  RVector<int> sampS_star;


  int q;
  int p;
  int K;
  int n;
  RVector<int> count;

  UpdateS_n(NumericVector S,NumericVector loglik, NumericVector X, NumericVector A,
          NumericVector beta,NumericVector u,
          IntegerVector sampS_star, int q, int p, int K, int n, IntegerVector count)
    : S(S), loglik(loglik), X(X), A(A), beta(beta),u(u),
      sampS_star(sampS_star), q(q), p(p), K(K), n(n), count(count){}


  void operator() (std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++){

      double S_star[q-1];
      double tf = 0;

      for(std::size_t l=0;l<(q-1);l++){
        S_star[l] = sampS_star[l+j*(q-1)];
        if( S_star[l] != S[l+(q-1)*j] ) {
          tf = 1;
        }
      }

      //compute log likelihood

      if(tf==1){
        double loglik_star = DLSS_compute_loglik_pll_n(X,A,beta,j+1,K,n,q,S_star);

        if(log(u[j]) < loglik_star-loglik[j]){
          for(std::size_t l=0;l<(q-1);l++){
            S[l+(q-1)*j] = S_star[l];
            count[j] = 1;
          }
          loglik[j] = loglik_star;
        }
      }
      else{
        count[j] = 1;
      }
    }
  }
};

struct UpdateA_n : public Worker {
  RVector<double> A;

  RVector<double> X;
  RVector<double> Y;

  double alpha_1;
  double alpha_0;

  int q;
  int p;
  int n;

  RVector<double> seed;
  RVector<double> seed2;

  UpdateA_n(NumericVector A, NumericVector X, NumericVector Y, double alpha_1,
          double alpha_0, int q, int p, int n, NumericVector seed, NumericVector seed2)
    : A(A), X(X), Y(Y), alpha_1(alpha_1), alpha_0(alpha_0),
      q(q), p(p), n(n), seed(seed), seed2(seed2) {}

  void operator() (std::size_t begin, std::size_t end) {

    for(std::size_t i = begin; i < end; i++){

      double b_i[q];
      double sum_b_i = 0;
      int ind = 0;
      double maxb = 0;

      for(std::size_t l=0; l<q; l++){
        b_i[l] = seed2[l];

      }

      for(std::size_t j=0;j<p;j++){
        ind = i+n*j;
        for(std::size_t l=0;l<q;l++){
          if(Y[ind]==l+1){
            b_i[l] += seed[ind];
            break;
          }
        }
      }

      for(int l=0; l<q; l++){
        if(b_i[l] > maxb){
          maxb = b_i[l];
        }
      }

      for(int l=0; l<q; l++){

        b_i[l] = exp( log(b_i[l]) - log(maxb) );

        sum_b_i += b_i[l];
      }

      for(std::size_t l=0;l<q;l++){
        A[i+n*l] = b_i[l] / sum_b_i;
      }

    }
  }
};



// [[Rcpp::export]]
void parallelDLSS_update_Y_n(NumericVector Y, NumericVector X, NumericVector A, NumericVector S,
                                    NumericVector seed, int q, int p, int n, int K){
  seed = Rcpp::runif(n*p,0,1);
  UpdateY_n updateY_n(Y, X, A, S, seed, q, p, n,K);
  parallelFor(0, p, updateY_n);

}


// [[Rcpp::export]]
void parallelDLSS_update_S_n(NumericVector S,NumericVector loglik, NumericVector X, NumericVector A, NumericVector beta,
                            NumericVector u, IntegerVector sampS_star, int q, int p, int K, int n, double lr, IntegerVector count){

  int sub = p*(q-1)*lr;

  if(sub ==0) {
    sub ++;
  }

  IntegerVector tag = Rcpp::sample(p*(q-1),sub,false) - 1;
  int loc = 0;

  for(int j = 0; j<p; j++){
    for(int l = 0; l<(q-1); l++){
      sampS_star[l + (q-1)*j] = (int)(S[l + (q-1)*j] + 1e-20);
    }
  }

  IntegerVector sam0 = Rcpp::sample(K,sub,true);

  for(int j = 0; j < sub; j++){
    loc = tag[j];
    sampS_star[loc] = sam0[j];
  }


  u = Rcpp::runif(p,0,1);

  UpdateS_n updateS_n(S, loglik, X, A, beta, u, sampS_star, q, p, K, n, count);
  parallelFor(0, p, updateS_n);

}


// [[Rcpp::export]]
void parallelDLSS_update_S_n_pickup(NumericVector S,NumericVector loglik, NumericVector X, NumericVector A, NumericVector beta,
                             NumericVector u, IntegerVector sampS_star, IntegerVector change, int q, int p, int K, int n, IntegerVector count){
  for(int j = 0; j<p; j++){
    if(change[j]==1){
      loglik[j] = -1e20;
    }

  }

  u = Rcpp::runif(p,0,1);

  UpdateS_n updateS_n(S, loglik, X, A, beta, u, sampS_star, q, p, K, n, count);
  parallelFor(0, p, updateS_n);

}

// [[Rcpp::export]]
void parallelDLSS_update_A_n(NumericVector A, NumericVector X, NumericVector Y, double alpha_1, double alpha_0,
                                    int q, int p, int n, NumericVector seed, NumericVector seed2){
  seed = -log( Rcpp::runif(n*p,0,1) + 1e-20 );
  for(int i=0; i<(q-1); i++){
    seed2[i] = R::rgamma(alpha_1+1e-20,1);
  }
  if(alpha_0==0){
    seed2[q-1] = 0;
  }
  else{
    seed2[q-1] = R::rgamma(alpha_0+1e-20,1);
  }

  UpdateA_n updateA_n(A, X, Y, alpha_1, alpha_0, q, p, n, seed,seed2);
  parallelFor(0, n, updateA_n);

}


using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List DLSS_gibbs_sampler_n(NumericVector X,NumericVector A0, NumericVector gamma, NumericVector alpha,
                                double lr, int total_iter = 1000, int burn_in = 100, int tune = 100, int show_step = 50) {
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A0.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];
  int K = (int) max(X);

  Rcpp::Dimension beta_dim(q-1,K);
  NumericVector beta(beta_dim);
  Rcpp::Environment MCMCpack("package:MCMCpack");
  Rcpp::Function rdirichlet = MCMCpack["rdirichlet"];

  Rcpp::Dimension S_dim(q-1,p);
  Rcpp::Dimension Y_dim(n,p);
  Rcpp::Dimension S_trace_dim(q-1,p,total_iter-burn_in);
  NumericVector S_trace(S_trace_dim);

  Rcpp::Dimension A_trace_dim(n,q,total_iter-burn_in);

  NumericVector A_trace(A_trace_dim);

  NumericVector S(S_dim);
  NumericVector Y(Y_dim);
  NumericVector A(A_dim);
  NumericVector loglik(p,-1e308);

  NumericVector seed(Y_dim);
  NumericVector seed2(q);

  NumericVector u(p);
  IntegerVector sampS_star(S_dim);
  IntegerVector count(p,0);
  int count0 = 0;

  for(int i=0;i<n;i++){
    for(int l=0;l<q;l++){
      A[i+A_dim[0]*l] = A0[i+A_dim[0]*l];
    }
  }


  for(int iter=0; iter< total_iter; iter++){

    for(int l=0; l<(q-1); l++){
      Rcpp::NumericVector beta_l = rdirichlet(Rcpp::_["n"] = K,Rcpp::_["alpha"] = gamma);
      for(int k=0; k<K; k++){
        beta[l*K+k] = beta_l[k];
      }
    }

    //update S
    parallelDLSS_update_S_n(S, loglik, X, A, beta,u, sampS_star, q, p, K, n,lr, count);

    count0 = 0;
    for(int kk = 0; kk < p; kk++){
      count0 = count0 + count[kk];
      count[kk] = 0;
    }

    //std::cout << lr << " " << 1.0*count0/p << std::endl;

    //update Y
    parallelDLSS_update_Y_n(Y, X, A, S, seed, q, p, n,K);

    //update A
    parallelDLSS_update_A_n(A, X, Y, alpha[0],alpha[1], q, p, n, seed, seed2);

    if(iter >= burn_in){
      //copy A_trace
      for(int i=0;i<n;i++){
        for(int l=0; l<q; l++){
          A_trace[i+A_dim[0]*l+(iter-burn_in)*A_dim[1]*A_dim[0]] = A[i+A_dim[0]*l];
        }
      }
      //copy S_trace
      for(int l=0;l<(q-1);l++){
        for(int j=0; j<p; j++){
          S_trace[l+S_dim[0]*j+(iter-burn_in)*S_dim[1]*S_dim[0]] = S[l+S_dim[0]*j];
        }
      }
    }

    if(iter < tune){
      if( 1.0*count0/p < 0.1 ){
        lr = lr * 0.9;
      }
      else if( 1.0*count0/p > 0.2){
        lr = lr * 1.1;
        if(lr > 1) {
          lr = 1;
        }
      }
    }


    if(iter%show_step==0){
      std::cout << "iter " << iter << std::endl;
    }
  }

  return Rcpp::List::create(Named("X")=X, Named("A")=A_trace,
                            Named("S") = S_trace,
                            Named("K")=K,
                            Named("logLik") = loglik,
                            Named("rate") = lr
                            );
}


// [[Rcpp::export]]
Rcpp::List DLSS_gibbs_sampler_n_pickup(NumericVector X,NumericVector A0,NumericVector loglik, NumericVector S0, NumericVector gamma, NumericVector alpha,
                                double noise, NumericVector S1, IntegerVector change, double lr, int total_iter = 1000, int burn_in = 100, int tune = 100, int show_step = 50) {
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A0.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];
  int K = (int) max(X);

  Rcpp::Dimension beta_dim(q-1,K);
  NumericVector beta(beta_dim);
  Rcpp::Environment MCMCpack("package:MCMCpack");
  Rcpp::Function rdirichlet = MCMCpack["rdirichlet"];

  Rcpp::Dimension S_dim(q-1,p);
  Rcpp::Dimension Y_dim(n,p);
  Rcpp::Dimension S_trace_dim(q-1,p,total_iter-burn_in);
  NumericVector S_trace(S_trace_dim);

  Rcpp::Dimension A_trace_dim(n,q,total_iter-burn_in);

  NumericVector A_trace(A_trace_dim);

  NumericVector S(S_dim);
  NumericVector Y(Y_dim);
  NumericVector A(A_dim);

  NumericVector seed(Y_dim);
  NumericVector seed2(q);

  NumericVector u(p);
  IntegerVector sampS_star(S_dim);

  IntegerVector count(S_dim,0);
  int count0 = 0;

  for(int i=0;i<n;i++){
    for(int l=0;l<q;l++){
      A[i+A_dim[0]*l] = A0[i+A_dim[0]*l];
    }
  }

  for(int l=0; l <(q-1); l++){
    for(int j=0; j<p; j++){
      S[l+S_dim[0]*j] = S0[l+S_dim[0]*j];
      sampS_star[l+S_dim[0]*j] = (int)( S1[l+S_dim[0]*j] + 1e-20);
    }
  }

  if(noise>0){
    for(int l=0; l<(q-1); l++){
      Rcpp::NumericVector beta_l = rdirichlet(Rcpp::_["n"] = K,Rcpp::_["alpha"] = gamma);
      for(int k=0; k<K; k++){
        beta[l*K+k] = beta_l[k];
      }
    }
    parallelDLSS_update_S_n_pickup(S, loglik, X, A, beta,u, sampS_star, change, q, p, K, n, count);
  }

  for(int iter=0; iter< total_iter; iter++){

    for(int l=0; l<(q-1); l++){
      Rcpp::NumericVector beta_l = rdirichlet(Rcpp::_["n"] = K,Rcpp::_["alpha"] = gamma);
      for(int k=0; k<K; k++){
        beta[l*K+k] = beta_l[k];
      }
    }

    //update S
    parallelDLSS_update_S_n(S, loglik, X, A, beta,u, sampS_star, q, p, K, n,lr,count);
    count0 = 0;
    for(int kk = 0; kk < p; kk++){
      count0 = count0 + count[kk];
      count[kk] = 0;
    }

    std::cout << lr << " " << 1.0*count0/p << std::endl;


    //update Y
    parallelDLSS_update_Y_n(Y, X, A, S, seed, q, p, n,K);

    //update A
    parallelDLSS_update_A_n(A, X, Y, alpha[0],alpha[1], q, p, n, seed, seed2);

    if(iter >= burn_in){
      //copy A_trace
      for(int i=0;i<n;i++){
        for(int l=0; l<q; l++){
          A_trace[i+A_dim[0]*l+(iter-burn_in)*A_dim[1]*A_dim[0]] = A[i+A_dim[0]*l];
        }
      }
      //copy S_trace
      for(int l=0;l<(q-1);l++){
        for(int j=0; j<p; j++){
          S_trace[l+S_dim[0]*j+(iter-burn_in)*S_dim[1]*S_dim[0]] = S[l+S_dim[0]*j];
        }
      }
    }
    if(iter < tune){
      if( 1.0*count0/p < 0.1 ){
        lr = lr * 0.9;
      }
      else if( 1.0*count0/p > 0.2){
        lr = lr * 1.1;
        if(lr > 1) {
          lr = 1;
        }
      }
    }


    if(iter%show_step==0){
      std::cout << "iter " << iter << std::endl;
    }
  }

  return Rcpp::List::create(Named("X")=X, Named("A")=A_trace,
                            Named("S") = S_trace,
                            Named("K")=K,
                            Named("logLik") = loglik,
                            Named("rate") = lr
                            );
}


