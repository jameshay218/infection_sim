#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double single_strain_model(double mu, double tp, double m, double y0, double ti, double t){
  double y = y0;
  if(ti == -1){
    return(y);
  }
  else {
    
    if(t <= ti){
      y = 0;
    }
    else if(t > ti && t <= (ti+tp)){
      y = (mu/tp)*t - (mu/tp)*ti;
    }
    else {
      y = -m*t + m*(ti+tp) + mu;
    }
    y = y + y0;
  }
  return(y);
}


//[[Rcpp::export]]
NumericMatrix single_infection_model(NumericVector mu_pars, NumericVector tp_pars, NumericVector m_pars, NumericVector y0_pars, int ti, NumericVector t){
  NumericMatrix y(t.size(), mu_pars.size());

  for(int j = 0;j < t.size();++j){
    for(int i =0; i < mu_pars.size();++i){
      y(j,i) = single_strain_model(mu_pars(i),tp_pars(i),m_pars(i),y0_pars(i),ti,t(j));
    }
  }

  return(y);
}

//[[Rcpp::export]]
NumericMatrix matrix_add(NumericMatrix a, NumericMatrix b){
  NumericMatrix c(a.nrow(),a.ncol());
  for(int i=0;i < c.nrow();++i){
    for(int j = 0; j < c.ncol();++j){
      c(i,j) = a(i,j) + b(i,j);
    }
  }
  return(c);
}


//[[Rcpp::export]]
NumericMatrix array_to_matrix(NumericVector a){
  NumericMatrix b(a.size()/3,a.size()/3);

  b(0,0) = a(0);
  b(0,1) = a(1);
  b(0,2) = a(2);
  b(1,0) = a(3);
  b(1,1) = a(4);
  b(1,2) = a(5);
  b(2,0) = a(6);
  b(2,1) = a(7);
  b(2,2) = a(8);

  return(b);
}


//[[Rcpp::export]]
 NumericMatrix multiple_infection_model(NumericMatrix mu_matrix, NumericMatrix tp_matrix, NumericMatrix m_matrix, NumericVector y0_pars, NumericVector ti_pars, NumericVector t){
   NumericMatrix y(t.size(),ti_pars.size());
   for(int i = 0;i<ti_pars.size();++i){
     y = matrix_add(y,single_infection_model(mu_matrix(i,_),tp_matrix(i,_),m_matrix(i,_),y0_pars, ti_pars(i),t));
   }

   return(y);
 }



//[[Rcpp::export]]
void print_matrix(NumericMatrix a){
  for(int i = 0; i < a.nrow();++i){
    for(int j = 0;j < a.ncol();++j){
      std::cout << a(i,j) << " ";
    }
    std:: cout << std::endl;
  }
}



//[[Rcpp::export]]
double optim_function(NumericVector params, NumericMatrix infection_times, NumericMatrix t0_titres, NumericMatrix t1_titres, double final_time){
  double sd = 1.0;
  int N = t0_titres.nrow();
  double ln = 0;
  int index = 0;

  NumericVector mu_pars(9);
  NumericVector tp_pars(9);
  NumericVector m_pars(9);
  NumericVector times(1);
  NumericVector single_tmp1(1);
  NumericVector lns;

  NumericMatrix mu_matrix(3,3);
  NumericMatrix tp_matrix(3,3);
  NumericMatrix m_matrix(3,3);
  NumericMatrix tmp(1,3);
  NumericMatrix titres(N,3);

  times[0] = final_time;

  for(int i =0;i <9;++i){
    mu_pars[index++] = params[i];
  }
  mu_matrix = array_to_matrix(mu_pars);

  // index = 0;
  //for(int i =9;i <18;++i){
   // tp_pars[index++] = params[i];
  //}
  //tp_matrix = array_to_matrix(tp_pars);

  index = 0;
  for(int i =9;i <18;++i){
    m_pars[index++] = params[i];
  }
  m_matrix = array_to_matrix(m_pars);  

  tp_pars = NumericVector::create(21,21,21,21,21,21,21,21,21);
  tp_matrix = array_to_matrix(tp_pars);

  for(int i=0;i<N;++i){
    tmp = multiple_infection_model(mu_matrix, tp_matrix, m_matrix, t0_titres(i,_),infection_times(i,_),times);
    titres(i,_) = tmp(0,_);
  }
    
  for(int i =0;i < titres.nrow();++i){
    for(int j = 0;j < titres.ncol();++j){
      single_tmp1[0] = titres(i,j);
      lns = dnorm(single_tmp1,t1_titres(i,j),sd,true);
      ln = ln +  std::accumulate(lns.begin(),lns.end(),0.0);
    }
  }

  return(-ln);  
  
}

