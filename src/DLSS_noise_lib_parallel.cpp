#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <chrono>
#include <ctime>
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
  int n;
  int K;

  UpdateY_n(NumericVector Y, NumericVector X, NumericVector A, NumericVector S,
            NumericVector seed, int q, int n, int K)
    : Y(Y), X(X), A(A), S(S), seed(seed),
      q(q), n(n), K(K){}


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

        tmp = seed[i+n*j] * sum_prob_Y;
        acm = 0;

        for(std::size_t l = 0; l<q; l++){
          acm += prob_Y[l];
          if(acm >= tmp){
            Y[ind] = l + 1;
            break;
          }
        }
      }
    }
  }
};

inline double DLSS_compute_loglik_pll_n(RVector<double> X,RVector<double> A,
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

  return loglik;
}


inline double s_density_n(RVector<double> X, RVector<double> A,
                          int j, int K, int n, int q,
                          double *S_j){

  j = j - 1;
  double loglik = 0.0;
  double temp = 0.0;

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

  return loglik;
}



struct UpdateS_n : public Worker {

  RVector<double> S;

  RVector<double> X;
  RVector<double> A;

  RVector<double> beta;
  RVector<int> group;

  int q;
  int p;
  int K;
  int G;
  int n;
  int l0;

  RVector<double> logprob;
  RVector<double> tmp;
  RVector<double> seed;

  UpdateS_n(NumericVector S,NumericVector X, NumericVector A, NumericVector beta, IntegerVector group,
            int q, int p, int K, int G, int n, int l0, NumericVector logprob, NumericVector tmp,
            NumericVector seed)
    : S(S), X(X), A(A), beta(beta), group(group), q(q), p(p), K(K), G(G), n(n),
      l0(l0), logprob(logprob),tmp(tmp), seed(seed){}

  void operator() (std::size_t begin, std::size_t end) {

    for(std::size_t j0 = begin; j0 < end; j0++){

      double tmp0[K];
      double Xind = 0;
      double tmp1 = 0;
      double logprob2[K];
      double prob[K];
      double maxprob;
      double oldS;
      double sum_prob;
      double acm;
      double tmp2;
      int g0 = 0;

      g0 = group[j0];

      for(int k = 0; k < K; k++){
        logprob[k + j0 * K] = log(beta[l0+g0*(q-1)+k*(q-1)*G] + 1e-20);
      }

      for(int i = 0; i < n; i++){

        Xind = X[i+n*j0];
        if(Xind==S[l0+(q-1)*j0]){
          tmp1 = tmp[i+n*j0] - A[i+n*l0] + 1e-10;
        }
        else{
          tmp1 = tmp[i+n*j0];
        }

        for(int k = 0; k < K; k++){
          tmp0[k] = tmp1;
        }

        tmp0[(int)(Xind+1e-20) - 1] += A[i+n*l0];

        for(int k = 0; k < K; k++){
          logprob[k + j0 * K] += log(tmp0[k]);
        }
      }

      for(int k = 0; k < K; k++){
        logprob2[k] = logprob[k + j0* K];
      }

      maxprob = -1e300;
      for(int k = 0; k < K; k++){
        if(logprob2[k] > maxprob){
          maxprob = logprob2[k];
        }
      }

      sum_prob = 0;
      for(int k = 0; k < K; k++){
        logprob2[k] = logprob2[k] - maxprob;
        prob[k] = exp(logprob2[k]);
        sum_prob += prob[k];
      }

      oldS = S[l0+(q-1)*j0];

      tmp2 = seed[l0+(q-1)*j0] * sum_prob;
      acm = 0;

      for(int k = 0; k < K; k++){
        acm += prob[k];
        if(acm > tmp2){
          S[l0+(q-1)*j0] = k + 1;
          break;
        }
      }

      if(oldS!=S[l0+(q-1)*j0]){

        for(int i = 0; i < n; i++){

          if(X[i + n*j0]==oldS){
            tmp[i+n*j0] -= A[i+n*l0];
          }
          else if(X[i + n*j0]==S[l0+(q-1)*j0]){
            tmp[i+n*j0] += A[i+n*l0];
          }
        }
      }

    }

  }
};


struct UpdateS_n_z : public Worker {

  RVector<double> S;

  RVector<double> X;
  RVector<double> Y;

  RVector<double> beta;
  RVector<int> group;

  int q;
  int p;
  int K;
  int G;
  int n;

  RVector<double> logprob;
  RVector<double> seed;

  UpdateS_n_z(NumericVector S,NumericVector X, NumericVector Y, NumericVector beta, IntegerVector group,
            int q, int p, int K, int G, int n, NumericVector logprob, NumericVector seed)
    : S(S), X(X), Y(Y), beta(beta), group(group), q(q), p(p), K(K), G(G), n(n),
      logprob(logprob), seed(seed){}

  void operator() (std::size_t begin, std::size_t end) {

    for(std::size_t j0 = begin; j0 < end; j0++){

      double tmp0[K];
      double tag[q];
      double tag0 = 0;
      double prob[K];
      double maxprob;
      double sum_prob;
      double acm;
      double tmp2;
      int g0 = 0;
      int ind = 0;
      int l0 = 0;

      g0 = group[j0];

      for(int l = 0; l < q; l++){
        tag[l] = 0;
      }

      for(int i = 0; i < n; i++){

        ind = i + j0 * n;
        l0 = Y[ind] - 1;

        if(l0 < (q-1) ){

          if(tag[l0]==0){
            S[l0+(q-1)*j0] = X[ind];
            tag[l0] = 1;
            tag0++;
          }

          if(tag0==(q-1) ){
            break;
          }
        }
      }

      if(tag0!=(q-1)){
        for(int l0 = 0; l0 < q; l0++){
          if(tag[l0]==0){
            for(int k = 0; k < K; k++){
              logprob[k + j0 * K] = log(beta[l0+g0*(q-1)+k*(q-1)*G] + 1e-200);
            }

            maxprob = -1e300;
            for(int k = 0; k < K; k++){
              if(logprob[k] > maxprob){
                maxprob = logprob[k];
              }
            }

            sum_prob = 0;
            for(int k = 0; k < K; k++){
              logprob[k] = logprob[k] - maxprob;
              prob[k] = exp(logprob[k]);
              sum_prob += prob[k];
            }

            tmp2 = seed[l0+(q-1)*j0] * sum_prob;
            acm = 0;

            for(int k = 0; k < K; k++){
              acm += prob[k];
              if(acm > tmp2){
                S[l0+(q-1)*j0] = k + 1;
                break;
              }
            }

          }

        }

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
      q(q), p(p), n(n), seed(seed), seed2(seed2){}

  void operator() (std::size_t begin, std::size_t end) {

    for(std::size_t i = begin; i < end; i++){

      double b_i[q];
      double sum_b_i = 0;
      int ind = 0;
      double maxb = 1e-200;

      for(std::size_t l=0; l<q; l++){
        b_i[l] = seed2[l+i*q];
      }

      for(std::size_t j=0;j<p;j++){

        ind = i+n*j;

        b_i[(int)(Y[ind]+1e-20)-1] += seed[i+n*j];

      }
/*
      for(int l=0; l<q; l++){
        if(b_i[l] > maxb){
          maxb = b_i[l];
        }
      }
*/
      for(int l=0; l<q; l++){

      //  b_i[l] = exp( log(b_i[l]) - log(maxb) );

        sum_b_i += b_i[l];
      }

      for(std::size_t l=0;l<q;l++){
        A[i+n*l] = b_i[l] / sum_b_i;
      }

    }
  }
};

struct Cal_tmp : public Worker {

  RVector<double> S;

  RVector<double> X;
  RVector<double> A;

  int q;
  int p;
  int n;
  int K;

  RVector<double> tmp;

  Cal_tmp(NumericVector S,NumericVector X, NumericVector A,
          int q, int p, int n, int K, NumericVector tmp)
    : S(S), X(X), A(A),q(q), p(p), n(n), K(K),tmp(tmp){}

  void operator() (std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++){

      int ind = 0;

      for(int i = 0; i < n; i++){

        ind = i + n*j;
        tmp[ind] = A[i+n*(q-1)]/K;

        for(int l = 0; l < (q-1); l++){
          if(X[ind]==S[l+(q-1)*j]){
            tmp[ind] += A[i+n*l];
          }
        }

      }

    }

  }
};

// [[Rcpp::export]]
void parallelDLSS_update_Y_n(NumericVector Y, NumericVector X, NumericVector A, NumericVector S,
                             NumericVector seed,int q, int p, int n, int K){
  seed = Rcpp::runif(n*p,0,1);

  UpdateY_n updateY_n(Y, X, A, S, seed, q, n,K);
  parallelFor(0, p, updateY_n);

}

// [[Rcpp::export]]
void parallelDLSS_update_Y_n_z(NumericVector Y, NumericVector X, NumericVector S,
                               NumericVector seed,int q, int p, int n, int K, double alpha_1, double alpha_0){
  seed = Rcpp::runif(n*p,0,1);

  NumericMatrix countY(n,q);
  NumericVector logprob(q);
  NumericVector prob(q);
  double logga = 0;
  double maxprob = 0;
  double sum_prob = 0;
  double tmp2 = 0;
  double acm = 0;
  NumericVector logga1(q);
  NumericVector candidate(q);

  for(int i = 0; i < n; i++){
    for(int l = 0; l < q; l++){
      countY(i,l) = 0;
    }
  }

  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      countY(i,Y[i+j*n]-1)++;
    }
  }

  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){

      countY(i,Y[i+j*n]-1)--;

      for(int l = 0; l < (q-1); l++){
        logga1(l) = lgamma(countY(i,l)+alpha_0);
        candidate(l) = 0;
      }
      logga1(q-1) = lgamma(countY(i,q-1)+alpha_1);
      candidate(q-1) = 0;

      for(int l = 0; l < q; l++){
        logga = 0;
        for(int l2 = 0; l2 < q; l2++){
          if(l2==l){
            logga += lgamma(countY(i,l2)+1);
          }else{
            logga += logga1(l2);
          }

        }

        if(l==(q-1)){
          logprob(q-1) = log(1.0/K) - logga;
          candidate(l) = 1;
        }else if(X[i+j*n]==S[l+(q-1)*j]){
          logprob(l) = -logga;
          candidate(l) = 1;
        }
      }

      maxprob = -1e300;
      for(int l = 0; l < q; l++){
        if(candidate(l) == 1){
          if(logprob(l) > maxprob){
            maxprob = logprob(l);
          }
        }
      }


      sum_prob = 0;
      for(int l = 0; l < q; l++){
        if(candidate(l) == 1){
          logprob(l) = logprob(l) - maxprob;
          prob(l) = exp(logprob(l));
          sum_prob += prob(l);
        }
      }

      tmp2 = seed[i+n*j] * sum_prob;

      acm = 0;
      for(int l = 0; l < q; l++){
        if(candidate(l) == 1){
          acm += prob(l);
          if(acm > tmp2){
            Y[i+n*j] = l + 1;
            break;
          }
        }
      }

      countY(i,Y[i+j*n]-1)++;

    }
  }


}



// [[Rcpp::export]]
void update_S_n(NumericVector S,NumericVector X, NumericVector A, NumericVector beta,
                IntegerVector group, int q, int p, int K, int G, int n){

  NumericVector seed = Rcpp::runif((q-1)*p,0,1);

  NumericVector logprob(K*p,0.0);
  NumericVector tmp(n*p);

  Cal_tmp cal_tmp(S, X, A, q, p, n, K,tmp);
  parallelFor(0, p, cal_tmp);

  for(int l = 0; l < (q-1); l++){

    UpdateS_n updateS_n(S, X, A,beta, group, q, p, K, G,  n, l, logprob,tmp,seed);
    parallelFor(0, p, updateS_n);

  }
}


// [[Rcpp::export]]
void update_S_n_z(NumericVector S,NumericVector X, NumericVector Y, NumericVector beta,
                IntegerVector group, int q, int p, int K, int G, int n){

  NumericVector seed = Rcpp::runif((q-1)*p,0,1);

  NumericVector logprob(K*p,0.0);

  UpdateS_n_z updateS_n_z(S, X, Y, beta, group, q, p, K, G,  n, logprob,seed);
  parallelFor(0, p, updateS_n_z);

}


// [[Rcpp::export]]
void parallelDLSS_update_A_n(NumericVector A, NumericVector X, NumericVector Y, double alpha_1, double alpha_0,
                                    int q,  int p, int n, NumericVector seed, NumericVector seed2){

  seed = -log( Rcpp::runif(n*p,0,1) + 1e-200 );

  for(int i=0; i < n; i++){
    for(int l=0; l<(q-1); l++){
      seed2[l+i*q] = R::rgamma(alpha_1+1e-20,1);
    }
    if(alpha_0==0){
      seed2[q-1+i*q] = 0;
    }
    else{
      seed2[q-1+i*q] = R::rgamma(alpha_0+1e-20,1);
    }
  }

  UpdateA_n updateA_n(A, X, Y, alpha_1, alpha_0, q, p, n, seed,seed2);
  parallelFor(0, n, updateA_n);

}


// [[Rcpp::export]]
void update_beta_n(NumericVector beta, NumericVector S, IntegerVector group,
                        double beta_0, int q, int p, int K, int G){

  NumericVector seed(q*p);
  NumericVector seed2(K);
  int g0 = 0;

  seed = -log( Rcpp::runif(q*p,0,1) + 1e-200 );

  NumericVector b_l(K*G);
  double sum_b_l = 0;
  int ind = 0;
  double maxb = 0;

  for(int l=0; l < (q-1); l++){

    for(int g=0; g < G; g++){

      for(int k=0; k<K; k++){
        seed2[k] = R::rgamma(beta_0+1e-20,1);
      }

      for(int k=0; k<K; k++){
        b_l[k + g*K] = seed2[k];
      }
    }

    for(int j=0;j<p;j++){

      ind = l+(q-1)*j;

      g0 = group[j];

      b_l[ (int)(S[ind]+1e-20)-1 + g0*K ] += seed[ind];

    }

    for(int g=0; g < G; g++){

      sum_b_l = 0;
      for(int k=0; k<K; k++){
        sum_b_l += b_l[k + g*K];
      }

      for(int k=0; k<K; k++){
        beta[l+ g*(q-1)+ k*(q-1)*G] = b_l[k + g*K] / sum_b_l;
      }

    }
  }
}


using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List NLSS_gibbs_sampler_n(NumericVector X,NumericVector A0, NumericVector S0,NumericVector Y0,NumericVector beta0, IntegerVector group,
                                double gamma, NumericVector alpha, int kk,
                                int total_iter = 1000, int burn_in = 100, int thin = 10, int show_step = 50) {
  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A0.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];
  int K = (int) max(X);
  int G = (int)( max(group) ) + 1;

  int ind = 0;

  int tag = 0;

  int nchain = (total_iter-burn_in)/thin;

  Rcpp::Dimension beta_dim((q-1)*G,K);
  NumericVector beta(beta_dim);

  Rcpp::Environment MCMCpack("package:MCMCpack");
  Rcpp::Function rdirichlet = MCMCpack["rdirichlet"];

  Rcpp::Dimension S_dim(q-1,p);
  Rcpp::Dimension Y_dim(n,p);
  Rcpp::Dimension seed_dim(n,p);

  Rcpp::Dimension S_trace_dim(q-1,p,nchain);

  NumericVector S_trace(S_trace_dim);

  Rcpp::Dimension beta_trace_dim((q-1)*G,K,nchain);
  NumericVector beta_trace(beta_trace_dim);

  Rcpp::Dimension A_trace_dim(n,q,nchain);

  NumericVector A_trace(A_trace_dim);

  NumericVector S(S_dim);
  NumericVector Y(Y_dim);
  NumericVector A(A_dim);

  NumericVector seed(seed_dim);
  NumericVector seed2(q*n);

  NumericVector u(p);

  std::clock_t start;
  double duration;


  auto end = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  IntegerVector sampS_star(S_dim);

  for(int i=0;i<n;i++){
    for(int l=0;l<q;l++){
      A[i+A_dim[0]*l] = A0[i+A_dim[0]*l];
    }
  }

  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      Y[i+n*j] = Y0[i+j*n];
    }
  }

  for(int l=0; l <(q-1); l++){
    for(int j=0; j<p; j++){
      S[l+S_dim[0]*j] = S0[l+S_dim[0]*j];
    }
  }

  for(int l=0; l<(q-1); l++){
    for(int g=0; g<G; g++){
      for(int k=0; k<K; k++){
        ind = l+g*(q-1)+k*(q-1)*G;
        beta[ind] = beta0[ind];
      }
    }
  }

  for(int iter=1; iter<=total_iter; iter++){

    //start = std::clock();
/*
    if((iter-burn_in)%kk==0){
      update_S_n(S, X, A, beta, group, q, p, K, G, n);
    }else{
      update_S_n_z(S, X, Y, beta, group, q, p, K, G, n);
    }
*/

    update_S_n_z(S, X, Y, beta, group, q, p, K, G, n);
    parallelDLSS_update_Y_n(Y, X, A, S, seed, q, p, n,K);
    parallelDLSS_update_A_n(A, X, Y, alpha[0],alpha[1], q, p, n, seed, seed2);
    update_beta_n(beta,S, group, gamma, q, p, K,G);

    /*
    if(iter < burn_in){

      // parallelDLSS_update_Y_n_z(Y, X,  S, seed, q, p, n,K,alpha[0],alpha[1]);
    }
    else{

      update_S_n(S, X, A, beta, group, q, p, K, G, n);
      parallelDLSS_update_Y_n(Y, X, A, S, seed, q, p, n,K);

    }

*/
    if( (iter > burn_in)&&( (iter-burn_in)%thin==0 ) ){
      //copy A_trace
      for(int i=0;i<n;i++){
        for(int l=0; l<q; l++){
          A_trace[i+A_dim[0]*l+tag*A_dim[1]*A_dim[0]] = A[i+A_dim[0]*l];
        }
      }
      //copy S_trace
      for(int l=0;l<(q-1);l++){
        for(int j=0; j<p; j++){
          S_trace[l+S_dim[0]*j+tag*S_dim[1]*S_dim[0]] = S[l+S_dim[0]*j];
        }
      }

      //copy beta_trace
      for(int l=0;l<(q-1);l++){
        for(int g=0; g<G; g++){
          for(int k=0; k<K; k++){
             beta_trace[l+(q-1)*g+(q-1)*G*k+tag*(q-1)*G*K] = beta[l+(q-1)*g+(q-1)*G*k];
          }
        }
      }

      tag++;
    }

    if(iter%show_step==0){
      end = std::chrono::system_clock::now();
      end_time = std::chrono::system_clock::to_time_t(end);
      std::cout << "iter " << iter  << " " <<  std::ctime(&end_time) << std::endl;
    }
  }

  return Rcpp::List::create(Named("X")=X, Named("A")=A_trace,
                            Named("S") = S_trace,
                            Named("beta") = beta_trace,
                            Named("Y") = Y,
                            Named("K")=K
                            );
}


// [[Rcpp::export]]
NumericVector cal_beta_coef(NumericVector S, int K){
  Rcpp::Dimension out_dim = S.attr("dim");
  NumericVector out(out_dim);

  int q = out_dim[0];
  int p = out_dim[1];
  int n = out_dim[2];

  for(int i = 0; i < q; i++){
    for(int j = 0; j < p; j++){
      for(int k = 0; k < n; k++){


      }
    }
  }

}
