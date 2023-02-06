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

struct ModeS : public Worker {

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

  ModeS(NumericVector S,NumericVector X, NumericVector A, NumericVector beta, IntegerVector group,
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
      double tmp2;
      int g0 = 0;
      int k0 = 0;

      g0 = group[j0];
      for(int k = 0; k < K; k++){
        logprob[k + j0 * K] = log(beta[l0+g0*(q-1)+k*(q-1)*G] + 1e-20);
        //logprob[k + j0 * K] = 0;
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

      maxprob = logprob2[0];

      for(int k = 1; k < K; k++){
        if(logprob2[k] > maxprob){
          maxprob = logprob2[k];
          k0 = k;
        }
      }

      oldS = S[l0+(q-1)*j0];
      S[l0+(q-1)*j0] = k0 + 1;

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

struct UpdateS : public Worker {

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

  UpdateS(NumericVector S,NumericVector X, NumericVector A, NumericVector beta, IntegerVector group,
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
      double sum_prob[K];
      double tmp2;
      int g0 = 0;

      g0 = group[j0];

      for(int k = 0; k < K; k++){
        logprob[k + j0 * K] = log(beta[l0+g0*(q-1)+k*(q-1)*G] + 1e-20);
        //logprob[k+ j0*K] = 0;
      }

      for(int i = 0; i < n; i++){

        Xind = X[i+n*j0];
        if(Xind==S[l0+(q-1)*j0]){
          tmp1 = tmp[i+n*j0] - A[i+n*l0];
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

      prob[0] = exp(logprob2[0]- maxprob);
      sum_prob[0] = prob[0];
      for(int k = 1; k < K; k++){
        prob[k] = exp(logprob2[k]- maxprob);
        sum_prob[k] = sum_prob[k-1] + prob[k];
      }

      oldS = S[l0+(q-1)*j0];
      tmp2 = seed[l0+(q-1)*j0] * sum_prob[K-1];

      for(int k = 0; k < K; k++){
        if(sum_prob[k] > tmp2){
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

struct UpdateS_t2 : public Worker {

  RVector<double> S;
  RVector<double> X;
  RVector<double> A;
  RMatrix<double> sumA;

  int q;
  int p;
  int K;
  int n;

  RVector<double> seed;
  RMatrix<double> Sspace;
  RVector<double> states;
  RMatrix<double> stat_link;
  int ns;

  UpdateS_t2(NumericVector S,NumericVector X, NumericVector A, NumericMatrix sumA,
             int q, int p, int K, int n,
             NumericVector seed, NumericMatrix Sspace,
             NumericVector states,
             NumericMatrix stat_link, int ns)
    : S(S), X(X), A(A), sumA(sumA), q(q), p(p), K(K), n(n),
      seed(seed), Sspace(Sspace), states(states),
      stat_link(stat_link), ns(ns){}

  void operator() (std::size_t begin, std::size_t end) {

    for(std::size_t j0 = begin; j0 < end; j0++){

      double logprob[ns];
      double maxprob;
      double sum_prob[ns];
      int Xind[n];
      double tmp2;
      int ind;

      for(int s = 0; s < ns; s++){
        logprob[s] = 0;
      }

      for(int i = 0; i < n; i++){
        for(int j = 0; j < K; j++){
          if(states[j]==X[i+n*j0]){
            Xind[i] = j;
            break;
          }
        }
      }

      for(int i = 0; i < n; i++){
        for(int s = 0; s < ns; s++){
          ind = stat_link(s,Xind[i]);
          logprob[s] += sumA(i,ind-1);
        }
      }

      maxprob = -1e300;
      for(int s=0; s<ns; s++){
        if(logprob[s] > maxprob){
          maxprob = logprob[s];
        }
      }

      sum_prob[0] = exp(logprob[0] - maxprob);

      for(int s = 1; s < ns; s++){
        sum_prob[s] =  sum_prob[s-1] + exp(logprob[s] - maxprob);
      }

      tmp2 = seed[j0] * sum_prob[ns-1];

      for(int s = 0; s < ns; s++){
        if(sum_prob[s] > tmp2){
          for(int l=0; l < (q-1); l++){
            S[l+(q-1)*j0] = Sspace(s,l);
          }
          break;
        }
      }

    }
  }
};

struct UpdateY : public Worker {

  RVector<double> Y;

  RVector<double> X;
  RVector<double> A;
  RVector<double> S;
  RVector<double> seed;

  int q;
  int n;
  int K;

  UpdateY(NumericVector Y, NumericVector X, NumericVector A, NumericVector S,
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

struct UpdateA : public Worker {
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

  UpdateA(NumericVector A, NumericVector X, NumericVector Y, double alpha_1,
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

      for(int l=0; l<q; l++){

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
void update_S_seprt(NumericVector S, NumericVector X, NumericVector A, NumericVector beta,
                    IntegerVector group, int q, int p, int K, int G, int n){

  NumericVector seed = Rcpp::runif((q-1)*p,0,1);

  NumericVector logprob(K*p,0.0);
  NumericVector tmp(n*p);

  Cal_tmp cal_tmp(S, X, A, q, p, n, K,tmp);
  parallelFor(0, p, cal_tmp);

  for(int l = 0; l < (q-1); l++){
    UpdateS updateS(S, X, A,beta, group, q, p, K, G,  n, l, logprob,tmp,seed);
    parallelFor(0, p, updateS);
  }
}


void recursive_sum_subset2(NumericMatrix a, NumericMatrix sum_a,
                           std::vector<std::string>& eqn_names, int& num_k, int k,
                           int nonzero, NumericVector nz_c){

  Rcpp::Dimension a_dim = a.attr("dim");
  int n = a_dim[0];
  int q = a_dim[1];

  if(k==1){
    for(int i = 0; i< n; i++){
      sum_a(i,0) = a(i,0);
    }
    eqn_names[0] = "a" + std::to_string(k);
    num_k = 1;
    nz_c[0] = 1;
  }
  if(k>1){
    int num_k_1;
    int tag = 1;

    recursive_sum_subset2(a, sum_a, eqn_names, num_k_1 ,k-1,nonzero, nz_c);

    eqn_names[num_k_1] = "a" + std::to_string(k);
    for(int i = 0; i< n; i++){
      sum_a(i,num_k_1) = a(i,k-1);
    }
    nz_c[num_k_1] = 1;
    for(int l=0;l<num_k_1;l++){
      if(nz_c[l]<nonzero){
        for(int i = 0; i< n; i++){
          sum_a(i,num_k_1+tag) = sum_a(i,l) + a(i,k-1);
        }
        eqn_names[num_k_1+tag] = eqn_names[l] + "a" + std::to_string(k);
        nz_c[num_k_1+tag] = nz_c[l] + 1;
        tag++;
      }
    }
    num_k = num_k_1 + tag;
  }
}

// [[Rcpp::export]]
Rcpp::List sum_subset2(NumericMatrix a, int q, int n, int nonzero) {

  int num0 = 1;
  int num1 = 1;

  int num_q = 0;

  for(int i = 1; i <= nonzero; i++){
    num0 = 1;
    num1 = 1;

    for(int j = q; j>(q-i); j--){
      num0 *= j;
    }

    for(int j = i; j > 0; j--){
      num1 *= j;
    }

    num_q += num0/num1;
  }

  NumericMatrix sum_a(n,num_q);
  NumericMatrix sum_a0(n,1);
  NumericVector nz_c(num_q, 0.0);
  std::vector<std::string> eqn_names(num_q );
  std::vector<std::string> eqn_names0(num_q + 1);

  recursive_sum_subset2(a,sum_a,eqn_names,num_q,q,nonzero,nz_c);

  eqn_names0[0] = "a0";

  for(int i = 1; i < (num_q + 1); i++){
    eqn_names0[i] = eqn_names[i-1];
  }

  return Rcpp::List::create(Named("sum_a") = cbind(sum_a0,sum_a),
                            Named("eqns") = eqn_names0 );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sum_subset3(NumericMatrix a, int q, int n, int nonzero) {

  int num0 = 1;
  int num1 = 1;

  int num_q = 0;

  for(int i = 1; i <= nonzero; i++){
    num0 = 1;
    num1 = 1;

    for(int j = q; j>(q-i); j--){
      num0 *= j;
    }

    for(int j = i; j > 0; j--){
      num1 *= j;
    }

    num_q += num0/num1;
  }


  NumericMatrix sum_a(n,num_q);
  NumericMatrix sum_a0(n,1);
  NumericVector nz_c(num_q, 0.0);
  std::vector<std::string> eqn_names(num_q );
  std::vector<std::string> eqn_names0(num_q + 1);

  recursive_sum_subset2(a,sum_a,eqn_names,num_q,q,nonzero,nz_c);

  eqn_names0[0] = "a0";
  for(int i = 1; i < (num_q + 1); i++){
    eqn_names0[i] = eqn_names[i-1];
  }
  return cbind(sum_a0,sum_a);
}


// [[Rcpp::export]]
void update_S_joint(NumericVector S,NumericVector X, NumericMatrix A, int q, int p,
                    int K, int n,
                    NumericMatrix Sspace, NumericVector states,
                    NumericMatrix stat_link, int ns, int q0){

  NumericVector seed = Rcpp::runif(p,0,1);

  int num0 = 1;
  int num1 = 1;
  int num_a = 0;

  for(int i = 1; i <= q0; i++){
    num0 = 1;
    num1 = 1;

    for(int j = (q-1); j>(q-1-i); j--){
      num0 *= j;
    }

    for(int j = i; j > 0; j--){
      num1 *= j;
    }

    num_a += num0/num1;
  }

  num_a += 1;

  NumericMatrix sumA(n, 2*num_a);
  NumericMatrix sumA1(n, num_a);
  NumericVector sumAall(n, 0.0);

  NumericVector newS((q-1)*p);

  sumA1 = sum_subset3(A,q-1,n,q0);

  Rcpp::Dimension sumA1_dim = sumA1.attr("dim");

  for(int i = 0; i < n; i++){
    for(int j = 0; j < (q-1); j++){
      sumAall(i) = sumAall(i) + A(i,j);
    }
    for(int j = 0; j < num_a; j++){
      sumA(i,j) = log(sumA1(i,j) + A(i,q-1)/K);
    }
    for(int j = num_a; j < (2*num_a); j++){
      sumA(i,j) = log(sumAall(i)- sumA1(i,j-num_a) + A(i,q-1)/K);
    }
  }

  UpdateS_t2 updateS_t2(newS, X, A, sumA, q, p, K, n, seed, Sspace, states,
                        stat_link, ns);
  parallelFor(0, p, updateS_t2);

  for(int l=0; l < (q-1); l++){
    for(int j=0; j<p; j++){
      S[l+(q-1)*j] = newS[l+(q-1)*j];
    }
  }

}


// [[Rcpp::export]]
void findmode_S(NumericVector S,NumericVector X, NumericVector A, NumericVector beta,
                IntegerVector group, int q, int p, int K, int G, int n){

  NumericVector seed = Rcpp::runif((q-1)*p,0,1);

  NumericVector logprob(K*p,0.0);
  NumericVector tmp(n*p);

  Cal_tmp cal_tmp(S, X, A, q, p, n, K,tmp);
  parallelFor(0, p, cal_tmp);

  for(int l = 0; l < (q-1); l++){

    ModeS modeS(S, X, A,beta, group, q, p, K, G,  n, l, logprob,tmp,seed);
    parallelFor(0, p, modeS);

  }
}


// [[Rcpp::export]]
void update_A(NumericVector A, NumericVector X, NumericVector Y, double alpha_1, double alpha_0,
              int q, int p, int n, NumericVector seed, NumericVector seed2){

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

  UpdateA updateA(A, X, Y, alpha_1, alpha_0, q, p, n, seed,seed2);
  parallelFor(0, n, updateA);

}

// [[Rcpp::export]]
void update_Y(NumericVector Y, NumericVector X, NumericVector A, NumericVector S,
              NumericVector seed,int q, int p, int n, int K){
  seed = Rcpp::runif(n*p,0,1);

  UpdateY updateY(Y, X, A, S, seed, q, n,K);
  parallelFor(0, p, updateY);

}

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List NLSS_gibbs_sampler(NumericVector X,NumericVector A0, NumericVector S0,NumericVector Y0,
                              NumericMatrix stat_link, int joint,
                              NumericVector alpha, NumericMatrix Sspace,
                              NumericVector states,
                              int q0, int K,
                              int total_iter = 1000, int burn_in = 100,
                              int thin = 10, int show_step = 50) {

  Rcpp::Dimension X_dim = X.attr("dim");
  Rcpp::Dimension A_dim = A0.attr("dim");
  Rcpp::Dimension Sspace_dim = Sspace.attr("dim");

  int n = X_dim[0];
  int p = X_dim[1];
  int q = A_dim[1];
  int ns = Sspace_dim[0];

  int ind = 0;
  int tag = 0;

  int nchain = (total_iter-burn_in)/thin;

  Rcpp::Environment MCMCpack("package:MCMCpack");
  Rcpp::Function rdirichlet = MCMCpack["rdirichlet"];

  Rcpp::Dimension S_dim(q-1,p);
  Rcpp::Dimension Y_dim(n,p);
  Rcpp::Dimension seed_dim(n,p);

  Rcpp::Dimension S_trace_dim(q-1,p,nchain);

  NumericVector S_trace(S_trace_dim);
  NumericVector Y_trace(S_trace_dim);

  Rcpp::Dimension A_trace_dim(n,q,nchain);

  NumericVector A_trace(A_trace_dim);

  NumericVector S(S_dim);
  NumericVector Y(Y_dim);
  NumericMatrix A(A_dim);

  NumericVector seed(seed_dim);
  NumericVector seed2(q*n);

  NumericVector u(p);

  int Sint = 0;

  std::clock_t start;
  double duration;

  auto end = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  for(int i=0;i<n;i++){
    for(int l=0;l<q;l++){
      A(i,l) = A0[i+n*l];
    }
  }

  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      Y[i+n*j] = Y0[i+j*n];
    }
  }

  for(int l=0; l <(q-1); l++){
    for(int j=0; j<p; j++){
      S[l+S_dim[0]*j] = S0[l+(q-1)*j];
    }
  }


  if(joint){
    for(int iter=1; iter<=total_iter; iter++){

      //std::cout << "update S" << std::endl;
      update_S_joint(S, X, A, q, p, K,  n, Sspace, states, stat_link, ns, q0);

      update_Y(Y, X, A, S, seed, q, p, n,K);

      update_A(A, X, Y, alpha[0],alpha[1], q, p, n, seed, seed2);


      if( (iter > burn_in)&&( (iter-burn_in)%thin==0 ) ){
        //copy A_trace
        for(int i=0;i<n;i++){
          for(int l=0; l<q; l++){
            A_trace[i+A_dim[0]*l+tag*A_dim[1]*A_dim[0]] = A(i,l);
          }
        }
        //copy S_trace
        for(int l=0;l<(q-1);l++){
          for(int j=0; j<p; j++){
            Sint = (int)(S[l+S_dim[0]*j]+1e-20);
            S_trace[l+S_dim[0]*j+tag*S_dim[1]*S_dim[0]] = Sint;
            S_trace[l+S_dim[0]*j+tag*S_dim[1]*S_dim[0]] = S[l+S_dim[0]*j];
          }
        }

        for(int l=0;l<(q-1);l++){
          for(int j=0; j<p; j++){
            Y_trace[l+S_dim[0]*j+tag*S_dim[1]*S_dim[0]] = Y[l+S_dim[0]*j];
          }
        }

        tag++;
      }

      if(iter%show_step==0){
        end = std::chrono::system_clock::now();
        end_time = std::chrono::system_clock::to_time_t(end);
        std::cout << " iter " << iter  << " " <<  std::ctime(&end_time) << std::endl;

      }

    }
  }
/*
  else{
    for(int iter=1; iter<=total_iter; iter++){

      update_Y(Y, X, A, S, seed, q, p, n,K);
      update_A(A, X, Y, alpha[0],alpha[1], q, p, n, seed, seed2);
      update_beta(beta,S, group, gamma, q, p, K,G);

      if(iter==1){
        findmode_S(S, X, A, beta, group, q, p, K, G, n);
      }
      else{
        update_S_seprt(S, X, A, beta, group, q, p, K, G, n);
      }

      if( (iter > burn_in)&&( (iter-burn_in)%thin==0 ) ){
        //copy A_trace
        for(int i=0;i<n;i++){
          for(int l=0; l<q; l++){
            A_trace[i+A_dim[0]*l+tag*A_dim[1]*A_dim[0]] = A(i,l);
          }
        }
        //copy S_trace
        for(int l=0;l<(q-1);l++){
          for(int j=0; j<p; j++){
            S_trace[l+S_dim[0]*j+tag*S_dim[1]*S_dim[0]] = (int)(S[l+S_dim[0]*j]);
          }
        }

        for(int l=0;l<(q-1);l++){
          for(int j=0; j<p; j++){
            Y_trace[l+S_dim[0]*j+tag*S_dim[1]*S_dim[0]] = Y[l+S_dim[0]*j];
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
  }
*/
  return Rcpp::List::create(Named("X") = X, Named("A") = A_trace,
                            Named("S") = S_trace,
                            Named("Y") = Y_trace,
                            Named("K")= K);
}


