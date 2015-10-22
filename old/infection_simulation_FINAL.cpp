#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <Rmath.h>
using namespace std;
using namespace Rcpp;

int MAX_TITRE = 15;
double EPSILON = 0.01;


//[[Rcpp::export]]
double f(double k, double lambda){
  double final=0;
  if(k != MAX_TITRE) return(R::dpois(k,lambda,0));
  else{
    for(int x = MAX_TITRE; x < 50;++x){
      final += R::dpois(x, lambda,0);
    }
  }
  return(final);
}
//[[Rcpp::export]]
double likelihood_simple(double lambda, double obs){
  return((1-((MAX_TITRE+1)*EPSILON)/MAX_TITRE)*f(obs,lambda) + EPSILON/MAX_TITRE);
}

//[[Rcpp::export]]
double likelihood_pois(NumericMatrix test, NumericMatrix data){
  double ln = 0;
  for(int i = 0; i < test.nrow(); ++i){
    for(int j = 0; j < test.ncol(); ++j){
      if(test(i,j) < 0) test(i,j) = 0;
      double tmp = likelihood_simple(test(i,j),data(i,j));
      if(tmp < -100000) cout << test(i,j) << " " << data(i,j) << endl;
      ln += log(likelihood_simple(test(i,j), data(i, j)));
    }
  }
  return(ln);
}


//[[Rcpp::export]]
double likelihood_norm(NumericMatrix test, NumericMatrix data, double sd, bool log){
  double ln = 0;
  for(int i = 0; i < test.nrow(); ++i){
    for(int j = 0; j < test.ncol(); ++j){
      ln += R::dnorm(test(i,j), data(i, j), sd, true);
    }
  }
  return(ln);
}

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

  //ln = likelihood_norm(titres, t1_titres, sd, true);
  ln = likelihood_pois(titres, t1_titres);
  return(ln);
}


//[[Rcpp::export]]
double proposal_function(double current, double step){
  double update;
  double move;
  update = R::rnorm(current, 1);
  current = current + (current-update)*step;
  return(current);
   
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
NumericVector multi_proposalfunction(NumericVector current, NumericVector lower, NumericVector upper, NumericVector step){
  NumericVector x;
  double tmp, rv;
  x = clone(current);
  for(int i = 0; i < x.size(); ++i){
    tmp = toUnitScale(x[i], lower[i], upper[i]);
    rv = (R::runif(0,1)-0.5)*step[i];
    tmp = tmp + rv;
    if(tmp < 0) tmp = -tmp;
    if(tmp > 1) tmp = 2-tmp;
    x[i] = fromUnitScale(tmp, lower[i], upper[i]);
  }
  return(x);
} 


//[[Rcpp::export]]
double scaletuning(double step, double popt, double pcur){
  if(pcur>=1) pcur = 0.99;
  if(pcur<=0) pcur = 0.01;
  step = (step*R::qnorm(popt/2,0,1,1,0))/R::qnorm(pcur/2,0,1,1,0);
  //  if(step > 1) step = 1;
  return(step);
}

//[[Rcpp::export]]
string run_MCMC(
		NumericMatrix param_table,
		int iterations, 
		NumericMatrix data, 
		double popt, 
		int opt_freq, 
		int thin, 
		int burnin, 
		int adaptive_period, 
		string filename,
		int save_block,
		bool VERBOSE,
		NumericMatrix infection_times,
		NumericMatrix t0_titres,
		double final_time
		){

  ofstream csv_write;
  double TUNING_ERROR = 0.1;
  double newprobab, probab, difflike, single_weight,samp_weight;
  int no_recorded = 0;
  int sampno = 1;
  int j = 0;
  string mcmc_chain_file = filename + "_chain.csv";
  
  NumericVector startvalue(param_table.nrow());
  NumericVector step_sizes(param_table.nrow());
  NumericVector current_params(param_table.nrow());
  NumericVector sampnos(save_block);
  NumericVector lnlikes(save_block);
  NumericVector tempaccepted(param_table.nrow());
  NumericVector tempiter(param_table.nrow());
  NumericVector reset(param_table.nrow());
  NumericVector pcur(param_table.nrow());
  NumericVector proposal(param_table.nrow());
  NumericVector lower_bounds(param_table.nrow());
  NumericVector upper_bounds(param_table.nrow());
  NumericVector weights(param_table.nrow());

  NumericMatrix empty_chain(save_block, param_table.nrow());
  NumericMatrix chain(save_block, param_table.nrow());
  
  startvalue = param_table(_,1);
  lower_bounds = param_table(_,2);
  upper_bounds = param_table(_,3);
  step_sizes = param_table(_,4);

  // Convert given weights to scale 0-1 for sampling
  single_weight = 1/(accumulate(param_table(_,5).begin(), param_table(_,5).end(), 0.0));
  weights[0] = param_table(0,5) * single_weight;
  for(int f = 1; f < param_table.nrow();++f){
    weights[f] = param_table(f,5)*single_weight + weights[f-1];
  }

  // Open stream to csv file
  csv_write.open(mcmc_chain_file.c_str());
  
  // Get initial parameter values and likelihood
  current_params = clone(startvalue);
  probab = posterior(startvalue, data, infection_times, t0_titres, final_time);

  // Write first line to csv file
  csv_write << sampno << ",";
  for(NumericVector::iterator it = startvalue.begin();it != startvalue.end();++it){
    csv_write << *it << ",";
  }
  csv_write << probab << endl;
  sampno += 1;

  // Print outs if asked
  if(VERBOSE){
    if(opt_freq==0) cout << "Not running adaptive MCMC - opt_freq set to 0" << endl;
    else cout << "Adaptive MCMC -will adapt step size during burnin period" << endl;
  }
    
  // For each iteration, go through each parameter
  for(int ii= 0; ii < (iterations + adaptive_period + burnin); ++ii){
    // Get random parameter based on weighting
    j = 0;
    samp_weight = R::runif(0,1);
    while(samp_weight > weights[j]){
      j++;
    }
    // Propose new value for current parameter
    proposal = clone(current_params);
    //proposal(j) = proposalfunction(proposal(j),lower_bounds(j),upper_bounds(j),step_sizes(j));
    proposal(j) = proposal_function(proposal(j),step_sizes(j));
    
    // Calculate likelihood with this new set of parameters
    newprobab = posterior(proposal, data, infection_times, t0_titres, final_time);
    
    // Get difference between this and old likelihood
    difflike = newprobab - probab;
    
    // Metropolis acceptance
    if(R::runif(0,1) < exp(difflike) || difflike > 0){
      current_params = clone(proposal);
      probab = newprobab;
      tempaccepted(j) += 1;
    }
    tempiter(j) += 1; 
    
    // If an iteration to be saved
    if(sampno%thin ==0){
      // Write sampno, current params and probab to chain
      sampnos(no_recorded) = sampno;
      for(int x = 0; x != current_params.size(); ++x){
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
    if(opt_freq != 0 && ii < adaptive_period && (ii+1)%opt_freq==0){
      for(int x = 0; x < current_params.size();++x){
	pcur[x] = tempaccepted(x)/tempiter(x);
	if(pcur(x) < popt - (TUNING_ERROR*popt) || pcur(x) > popt + (TUNING_ERROR*popt)){
	  step_sizes(x) = scaletuning(step_sizes(x), popt, pcur(x));
	}
      }
      tempaccepted = clone(reset);
      tempiter = clone(reset);
    }
  }
  csv_write.close();
  return(mcmc_chain_file);
}
