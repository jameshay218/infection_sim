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

//[[Rcpp::export]]
NumericMatrix individual_sim(
			     NumericVector mu_pars, 
			     NumericVector tp_pars, 
			     NumericVector m_pars, 
			     NumericVector ti_pars, 
			     double  y0b,
			     double lower_titre_bound, 
			     NumericVector times
			     ){
  double y0 =  y0b;
  double final_t, t_i, mu, tp, m, tmp, max_t, tmp2;
  NumericVector::iterator t = times.begin();
  NumericVector Y(times.size());
  NumericMatrix out(times.size(),2);

  int j = 0;
  int no_infections = ti_pars.size();
  max_t = times[times.size()-1];
  
  for(int i = 1; i <=no_infections; ++i){
    if(i == no_infections){
      final_t = max_t;
    }
    else {
      final_t = ti_pars[i];
    }
    t_i = ti_pars[i-1];
    mu = mu_pars[i-1];
    tp = tp_pars[i-1];
    m = m_pars[i-1];
    
    while(t != times.end() && *t <= final_t){
      tmp = 0;
      if(*t <= t_i) tmp = y0;
      else if(*t > t_i && *t <= (t_i + tp)) tmp = (mu/tp)*(*t) - (mu/tp)*(t_i) + y0;
      else tmp = -m*((*t)) + m*(t_i + tp) + mu + y0;
      Y[j] = tmp;
      
      if(Y[j] < lower_titre_bound) Y[j] = lower_titre_bound;
      ++t;
      ++j;
    }
    
    if(i < no_infections){
      tmp2 = ti_pars[i];
      if(tmp2 <= t_i) y0 = y0;
      else if (tmp2 <= t_i + tp) y0 = (mu/tp)*tmp2 - (mu/tp)*t_i + y0;
      else y0 = m*(t_i + tp - tmp2) + mu + y0;
    }
   
    if(y0 < lower_titre_bound){
      y0 = lower_titre_bound;
    }
  }
  out(_,0) = times;
  out(_,1)=Y;
  return out;
}



//[[Rcpp::export]]
double obs_error(int actual, int obs, double S, double EA){
  int MAX_TITRE = 8;
  //  double S = 0.95;
  //  double EA = 0.02;
  if(actual == (MAX_TITRE) && obs == (MAX_TITRE)) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==0 && obs==0) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==obs) return(S);
  else if(actual == (obs + 1) || actual==(obs-1)) return(EA/2.0);
  return((1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
}


//[[Rcpp::export]]
double toUnitScale(double x, double min, double max){
  return((x-min)/(max-min));
}

//[[Rcpp::export]]
double fromUnitScale(double x, double min, double max){
  return(min + (max-min)*x);
}



//[[Rcpp::export]]
double proposalfunction(double current, double lower, double upper, double step){
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
double proposal_function(double current, double lower, double upper, double step){
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
double posterior(NumericVector ti_pars, NumericMatrix mu_pars, NumericMatrix tp_pars, NumericMatrix m_pars, NumericVector other_pars, NumericMatrix data,NumericVector time){
  // Sort infection times into correct order
  double y0 = other_pars[0];
  double S = other_pars[1];
  double EA = other_pars[2];
  double lower_bound = other_pars[3];
  NumericMatrix y(data.nrow(), 2);
  double ln = 0;

  NumericVector new_ti = clone(ti_pars);
  IntegerVector indices(ti_pars.size());
  NumericVector mus(ti_pars.size());
  NumericVector tps(ti_pars.size());
  NumericVector ms(ti_pars.size());
  new_ti = vector_sort(ti_pars);
  indices = vector_order(ti_pars) - 1;
  
  for(int i = 0; i < ti_pars.size(); ++i){
    mus = mu_pars(i,_);
    mus = mus[indices];
    tps = tp_pars(i,_);
    tps = tps[indices];
    ms = m_pars(i,_);
    ms = ms[indices];
    y = individual_sim(mus,tps,ms, new_ti,y0, lower_bound,time);
    for(int j = 0; j < y.nrow();++j){
      // cout << "y: " << y(j,1) << endl;
      //cout << "data: " << data(j,i) << endl;
      ln += log(obs_error(floor(y(j,1)),floor(data(j,i)),S,EA));
      //ln += R::dnorm(y(j,1),data(j,i),1,true);
      //ln += R::dpois(floor(data(j,1)),floor(y(j,i)),true);
    }
  }
  return(ln);
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


//[[Rcpp::export]]
double scaletuning(double step, double popt, double pcur){
  if(pcur>=1) pcur = 0.99;
  if(pcur<=0) pcur = 0.01;
  step = (step*R::qnorm(popt/2,0,1,1,0))/R::qnorm(pcur/2,0,1,1,0);
  if(step > 10) step = 9.99;
  if(step < 0) step = 0.01;
  return(step);
}


// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   bool success = false;
   arma::mat Y = arma::randn(n, ncols);
   arma::mat R;
   int j = 0;
   while(!success){
     success = arma::chol(R, sigma);
     if(!success){
       //cout << "Error in mvnorm" << endl;
       sigma += arma::eye(sigma.n_rows, sigma.n_rows)*1e-5;
       j++;
       cout << j << endl;
     }
   }
   return arma::repmat(mu, 1, n).t() + Y * R;
}


//[[Rcpp::export]]
arma::mat subset_mat(arma::mat a, arma::uvec b){
  return(a.submat(b,b));
}

//[[Rcpp::export]]
NumericVector mvnorm_proposal(NumericVector current, arma::mat covar_mat, NumericVector upper, NumericVector lower){
  NumericVector proposed(current.size());
  proposed = clone(current);
  NumericVector x(current.size());
  NumericVector zeroes(current.size());
  proposed = clone(current);
  x = mvrnormArma(1,zeroes,covar_mat);
  for(int i = 0; i < current.size();++i){
    proposed(i) = current(i) + x(i);
  }
  return(proposed);
}



//[[Rcpp::export]]
string run_MCMC(NumericVector startvalue,
		NumericMatrix data,
		NumericVector times,
		NumericMatrix param_table,
		int iterations, 
		double popt, 
		int opt_freq, 
		int thin, 
		int burnin, 
		int adaptive_period, 
		string filename,
		int save_block,
		NumericMatrix mu_pars,
		NumericMatrix tp_pars,
		NumericMatrix m_pars,
		NumericVector control_pars
		){
  ofstream csv_write; // Output stream for csv results
  double TUNING_ERROR = 0.1; // Allowable tuning error for adaptive step sizes
  double newprobab, probab, difflike; // Doubles for log-likelihoods and weights
  int no_recorded = 0; // Record number of recorded iterations
  int sampno = 1; // Record current sample number (note this may be different to no_recorded if thin != 1)
  int j = 0; // Record current parameter index
  string mcmc_chain_file = filename + "_chain.csv";

  // Vectors for parameters
  NumericVector step_sizes(startvalue.length());
  NumericVector current_params(startvalue.length());
  NumericVector sampnos(save_block);
  NumericVector lnlikes(save_block);
  NumericVector proposal(startvalue.length());
  vector<int> fixed;
  NumericVector log_proposal(startvalue.length());
  
  // Vectors for adaptive steps
  NumericVector tempaccepted(startvalue.length()); // Store total number of accepted proposals for each parameter
  NumericVector tempiter(startvalue.length()); // Store total number of proposals for each parameter
  NumericVector reset(startvalue.length()); // Vector of zeroes for reset of above vectors
  NumericVector pcur(startvalue.length()); // Vector of acceptance rates for each parameter

  /*double tempaccepted = 0;
    double tempiter = 0;
    double reset= 0;
    double pcur = 0;
  */

  NumericVector lower_bounds(param_table.nrow()); // Lower bounds for each parameter
  NumericVector upper_bounds(param_table.nrow()); // Upper bounds for each parameter

  // Matrices to store the parameter set at each time point
  NumericMatrix empty_chain(save_block, startvalue.length());
  NumericMatrix chain(save_block, startvalue.length());
  
  // Get control parameters from parameter table
  for(int i =0; i < startvalue.length();++i){
    if(param_table(i,1) == 0){
      fixed.push_back(i);
    }
    step_sizes[i] = param_table(i,4);
    lower_bounds[i] = param_table(i,2);
    upper_bounds[i] = param_table(i,3);
  }

  // Open stream to csv file
  csv_write.open(mcmc_chain_file.c_str());
  
  // Get initial parameter values and likelihood
  current_params = clone(startvalue);
  probab = posterior(current_params, mu_pars, tp_pars, m_pars, control_pars, data, times);


  // Write first line to csv file
  csv_write << sampno << ",";
  for(int i = 0; i < (startvalue.length()); ++i){
    csv_write << startvalue[i] << ",";
  }
  csv_write << probab << endl;
  sampno += 1;
  int jj = 0;


  //arma::mat cov_matrix(start_covmatrix.nrow(),start_covmatrix.ncol());
  //arma::mat cov_matrix_scaled(start_covmatrix.nrow(),start_covmatrix.ncol());

  //  cov_matrix = as<arma::mat>(start_covmatrix);
  //  cov_matrix_scaled = cov_matrix*scale;


  for(int i = 0; i < (iterations + adaptive_period + burnin); ++i){
    jj = randNumber(fixed.size());
    j = fixed[jj];
    proposal = clone(current_params);
    //proposal(j) = proposal_function(proposal(j),lower_bounds(j),upper_bounds(j),step_sizes(j));

    proposal(j) = proposal_function(proposal(j),lower_bounds(j), upper_bounds(j),step_sizes(j));
    //proposal = mvnorm_proposal(proposal, cov_matrix_scaled, upper_bounds, lower_bounds);
    newprobab = posterior(proposal, mu_pars, tp_pars, m_pars, control_pars, data, times);
    
    difflike=newprobab-probab;
    if(R::runif(0,1) < exp(difflike) || difflike > 0){
      current_params = clone(proposal);
      probab = newprobab;
      tempaccepted(j) += 1;
      //tempaccepted += 1;
    }
    tempiter(j) += 1; 
    //tempiter += 1;
      
    // If an iteration to be saved
    if(sampno%thin ==0){
      // Write sampno, current params and probab to chain
      sampnos(no_recorded) = sampno;
      for(int x = 0; x < (current_params.size()); ++x){
	chain(no_recorded,x) = current_params[x];
      }
      lnlikes(no_recorded) = probab;
      no_recorded++;
      
      if(no_recorded >= save_block){
	// Write chain to csv file
	for(int x = 0; x < chain.nrow(); ++x){
	  csv_write << sampnos(x) << ",";
	  for(int q =0; q < chain.ncol();++q){
	    csv_write << chain(x,q) << ",";
	  }
	  csv_write << lnlikes(x);
	  csv_write << endl;
	}
	chain = empty_chain;
	no_recorded = 0;
      }
    }
    sampno++;

    // If in the adaptive period and an update step
    if(opt_freq != 0 && i < adaptive_period && (i+1)%opt_freq==0){

      // Update step sizes for normal parameters
      for(int x = 0; x < fixed.size();++x){
	int jjj = fixed[x];
	// get current acceptance rate
	if(tempiter(jjj) > 0){
	  pcur[jjj] = tempaccepted(jjj)/tempiter(jjj);
	  // If outside acceptable range:
	  if(pcur(jjj) < popt - (TUNING_ERROR*popt) || pcur(jjj) > popt + (TUNING_ERROR*popt)){
	    step_sizes(jjj) = scaletuning(step_sizes(jjj), popt, pcur(jjj));
	    tempaccepted = clone(reset);
	    tempiter = clone(reset);
	    
	  }
	}
      }
      //      pcur = tempaccepted/tempiter;
      //scale = scaletuning(scale, popt, pcur);
      /*cout << scale << endl;
      arma::mat greb = cov(as<arma::mat>(chain));
      cov_matrix = (1-w)*cov_matrix + w*greb;
      cov_matrix_scaled = cov_matrix*scale;
      tempiter = tempaccepted = 0;*/
    }
  }
  csv_write.close();

  return(mcmc_chain_file);
}













//[[Rcpp::export]]
NumericMatrix OLD_individual_sim_OLD(NumericVector pars, NumericVector times){
  double lower_titre_bound = pars[0];
  double y0 = pars[1];
  double final_t, t_i, mu, tp, m, tmp, max_t, tmp2;
  NumericVector::iterator t = times.begin();
  NumericVector Y(times.size());
  NumericMatrix out(times.size(),2);

  int j = 0;
  int no_infections = (pars.size()-2)/4;
  max_t = times[times.size()-1];
  for(int i = 1; i <=no_infections; ++i){
    if(i == no_infections){
      final_t = max_t;
    }
    else {
      final_t = pars[4*i + 2];
    }

   t_i = pars[4*(i-1) + 2];
   mu = pars[4*(i-1) + 3];
   tp = pars[4*(i-1) + 4];
   m = pars[4*(i-1) + 5];

   while(t != times.end() && *t <= final_t){
     tmp = 0;
     if(*t <= t_i) tmp = y0;
     else if(*t > t_i && *t <= (t_i + tp)) tmp = (mu/tp)*(*t) - (mu/tp)*(t_i) + y0;
     else tmp = -m*((*t)) + m*(t_i + tp) + mu + y0;
     Y[j] = tmp;

     if(Y[j] < lower_titre_bound) Y[j] = lower_titre_bound;
     ++t;
     ++j;
   }

   if(i < no_infections){
     tmp2 = pars[4*i + 2];
     if(tmp2 <= t_i) y0 = y0;
     else if (tmp2 <= t_i + tp) y0 = (mu/tp)*tmp2 - (mu/tp)*t_i + y0;
     else y0 = m*(t_i + tp - tmp2) + mu + y0;
   }
   
   if(y0 < lower_titre_bound){
     y0 = lower_titre_bound;
   }
  }
out(_,0) = times;
  out(_,1)=Y;
  return out;
}
