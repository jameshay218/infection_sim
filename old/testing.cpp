#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <Rmath.h>
using namespace std;
using namespace Rcpp;


//[[Rcpp::export]]
double single_strain_model(double mu, double tp, double m, double ti, double t){
  if(ti == -1) return(0);
  if(t <= ti)  return(0);
  if(t > ti && t <= (ti+tp)) return((mu/tp)*t - (mu/tp)*ti);
  return(-m*t + m*(ti+tp) + mu);
}

//[[Rcpp::export]]
void multiple_infection_model(NumericMatrix mu_pars, NumericMatrix tp_pars, NumericMatrix m_pars, NumericMatrix ti_pars, double t,  NumericMatrix y){
  // For each individual
  for(int i =0; i < y.nrow();++i){
  // For each infection
    for(int j = 0;j < ti_pars.ncol();++j){
      // For each boost
      for(int x = 0; x < mu_pars.ncol();++x){
	y(i,x) += single_strain_model(mu_pars(j,x),tp_pars(j,x),m_pars(j,x),ti_pars(i,j),t);
      }
    }
  }
}


//[[Rcpp::export]]
double posterior(NumericVector params, NumericMatrix t1_titres, NumericMatrix infection_times, NumericMatrix t0_titres, double final_time){
  double sd = 1.0;
  int N = t0_titres.nrow();
  double ln = 0;
  NumericMatrix mu_matrix(3,3);
  NumericMatrix tp_matrix(3,3);
  NumericMatrix m_matrix(3,3);
  NumericMatrix titres(N,3);
  titres = clone(t0_titres);
  
  // Put mu params into matrix
  mu_matrix(0,0) = params[0];
  mu_matrix(0,1) = params[1];
  mu_matrix(0,2) = params[2];
  mu_matrix(1,0) = params[3];
  mu_matrix(1,1) = params[4];
  mu_matrix(1,2) = params[5];
  mu_matrix(2,0) = params[6];
  mu_matrix(2,1) = params[7];
  mu_matrix(2,2) = params[8];

  // Put tp params into matrix
  tp_matrix(0,0) = 21;
  tp_matrix(0,1) = 21;
  tp_matrix(0,2) = 21;
  tp_matrix(1,0) = 21;
  tp_matrix(1,1) = 21;
  tp_matrix(1,2) = 21;
  tp_matrix(2,0) = 21;
  tp_matrix(2,1) = 21;
  tp_matrix(2,2) = 21;

  // Put m params into matrix
  m_matrix(0,0) = params[9];
  m_matrix(0,1) = params[10];
  m_matrix(0,2) = params[11];
  m_matrix(1,0) = params[12];
  m_matrix(1,1) = params[13];
  m_matrix(1,2) = params[14];
  m_matrix(2,0) = params[15];
  m_matrix(2,1) = params[16];
  m_matrix(2,2) = params[17];

  multiple_infection_model(mu_matrix, tp_matrix, m_matrix, infection_times, final_time, titres);

  for(int i = 0; i < titres.nrow(); ++i){
    for(int j = 0; j < titres.ncol(); ++j){
      ln += R::dnorm(titres(i,j), t1_titres(i, j), sd, true);
    }
  }
  return(ln);
}

//[[Rcpp::export]]
NumericVector weight_test(NumericMatrix param_table){
  NumericVector weights(param_table.nrow());
  double single_weight;
  single_weight = 1/(accumulate(param_table(_,5).begin(), param_table(_,5).end(), 0.0));
  weights[0] = param_table(0,5) * single_weight;
  for(int i = 1; i < param_table.nrow();++i){
    cout << i << endl;
    weights[i] = param_table(i,5)*single_weight + weights[i-1];
  }
  return(weights);
}

//[[Rcpp::export]]
int sample_param(NumericVector weights){
  double samp_weight;
  int j = 0;
  samp_weight = R::runif(0,1);
  cout << samp_weight << endl;
  while(samp_weight > weights[j]){
    j++;
  }
  cout << weights[j] << endl;
  return(j);
}
//[[Rcpp::export]]
NumericVector test_copy(NumericMatrix param_table){
  NumericVector test(param_table.nrow());
  test = param_table(_,1);
  test[0] = 9999;
  return(test);
}

//[[Rcpp::export]]
double f1(double k, double lambda){
  return(R::dpois(k,lambda,true));
}
  
int MAX_TITRE = 15;
//[[Rcpp::export]]
double f(double k, double lambda){
  double final=0;
  if(k != MAX_TITRE) return(R::dpois(k,lambda,true));
  else{
    for(int x = MAX_TITRE; x < 50;++x){
      final += R::dpois(x, lambda,true);
    }
  }
  return(final);
}
