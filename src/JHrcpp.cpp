//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <iostream>
#include <fstream>
#include <Rmath.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace std;
using namespace Rcpp;

/* =============================================== */
/* --------------------------------- MATHS FUNCTION --------------------------------------------- */
/* =============================================== */
//[[Rcpp::export]]
NumericVector vector_sort(NumericVector x){
  NumericVector y = clone(x);
  return(y.sort());
}
//[[Rcpp::export]]
IntegerVector vector_order(NumericVector x){
  NumericVector y = clone(x).sort();
  return(match(y,x));
}

//[[Rcpp::export]]
double toUnitScale(double x, double min, double max){
  return((x-min)/(max-min));
}

//[[Rcpp::export]]
double fromUnitScale(double x, double min, double max){
  return(min + (max-min)*x);
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   bool success = false;
   arma::mat Y = arma::randn(n, ncols);
   arma::mat R;
   while(!success){
     success = arma::chol(R, sigma);
     if(!success){
       cout << "Error in mvnorm" << endl;
       sigma += arma::eye(sigma.n_rows, sigma.n_rows)*1e-3;
     }
   }
   return arma::repmat(mu, 1, n).t() + Y * R;
}

//[[Rcpp::export]]
arma::mat subset_mat(arma::mat a, arma::uvec b){
  return(a.submat(b,b));
}

// Random number from 0 to N
//[[Rcpp::export]]
int randNumber(const int n){return floor(unif_rand()*n);}


// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// it would appear that randomSample is still only a GNU g++ extension ?
#include <ext/algorithm>

// [[Rcpp::export]]
Rcpp::NumericVector randomSample(Rcpp::NumericVector a, int n) {
    // clone a into b to leave a alone
    Rcpp::NumericVector b(n);
    __gnu_cxx::random_sample(a.begin(), a.end(), 
                             b.begin(), b.end(), randWrapper);
    return b;
}



/* =============================================== */
/* --------------------------------- PROPOSAL FUNCTIONS -------------------------------------- */
/* =============================================== */
//[[Rcpp::export]]
double uvnorm_proposal(double current, double lower, double upper, double step){
  double update;
  double move;
  double new1;
  new1 = toUnitScale(current,lower,upper);
  
  do {
    new1 = toUnitScale(current,lower,upper);
    update = R::rnorm(0, 1);
    new1 = new1 + update*step;
  } while(new1 > 1 || new1 < 0);
  
  new1 = fromUnitScale(new1,lower,upper);
  
  return(new1);
}

//[[Rcpp::export]]
double uu_proposal(double current, double lower, double upper, double step){
  double x,rv,rtn;
  x = toUnitScale(current, lower, upper);
  rv = (R::runif(0,1)-0.5)*step;
  x = x + rv;
  if (x < 0) x = -x;
  if(x > 1) x = 2-x;

  rtn = fromUnitScale(x,lower,upper);
  return(rtn);
}


//[[Rcpp::export]]
NumericVector mvnorm_proposal(NumericVector current, arma::mat covar_mat, NumericVector upper, NumericVector lower, NumericVector fixed){
  NumericVector proposed(current.size());
  proposed = clone(current);
  NumericVector x(fixed.size());
  NumericVector zeroes(fixed.size());
  proposed = clone(current);
  int j = 0;
  int ii = 0;
  while(j==0){
    j = 1;
    x = mvrnormArma(1,zeroes,covar_mat);
    for(int i = 0; i < fixed.size();++i){
      ii = fixed[i];
      proposed(ii) = current(ii) + x(i);
    }  
  }
  return(proposed);
}


/* =============================================== */
/* ---------------------------------- AUXILIARY MCMC --------------------------------------------- */
/* =============================================== */

//[[Rcpp::export]]
double scaletuning(double step, double popt, double pcur){
  if(pcur>=1) pcur = 0.99;
  if(pcur<=0) pcur = 0.01;
  step = (step*R::qnorm(popt/2,0,1,1,0))/R::qnorm(pcur/2,0,1,1,0);
  if(step > 1) step = 0.99;
  if(step < 0) step = 0.01;
  return(step);
}

