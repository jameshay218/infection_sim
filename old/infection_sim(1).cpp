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

int MAX_TITRE = 15;
double EPSILON = 0.01;

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
double posterior1(NumericVector params, NumericMatrix t1_titres, NumericMatrix infection_times, NumericMatrix t0_titres, double final_t){
  double sd=1.0;
  int N = t0_titres.nrow();
  double ln = 0;
  
  NumericMatrix mu_matrix(1,3);
  NumericMatrix tp_matrix(1,3);
  NumericMatrix m_matrix(1,3);
  NumericMatrix titres(N,1);
  titres = clone(t0_titres);
  
  mu_matrix(0,0) = params[0];
  mu_matrix(0,1) = params[1];
  mu_matrix(0,2) = params[2];

  tp_matrix(0,0) = 21;
  tp_matrix(0,1) = 21;
  tp_matrix(0,2) = 21;
  
  m_matrix(0,0) = params[3];
  m_matrix(0,1) = params[4];
  m_matrix(0,2) = params[5];
    
 multiple_infection_model(mu_matrix, tp_matrix, m_matrix, infection_times, final_t, titres);

  ln = likelihood_norm(titres, t1_titres, sd, true);
  return(ln);
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
   cout << "Titre: " << titres(0,0) << endl;
   cout << "True titre: " << t1_titres(0,0) << endl;
  ln = likelihood_norm(titres, t1_titres, sd, true);
    cout << "ln: " << ln << endl;
  //ln = likelihood_pois(titres, t1_titres);
  return(ln);
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
double proposal_function(double current, double step){
  double update;
  double move;
  update = R::rnorm(current, 1);
  current = current + (current-update)*step;
  return(current);
}

//[[Rcpp::export]]
NumericMatrix infection_proposalfunction(NumericMatrix current, double lower, double upper, double step,int N){
  NumericMatrix x;
  double tmp, rv;
  NumericVector samples;
  samples = runif(N)*current.nrow();
  x = clone(current);
  int index = 0;
  for(int i = 0; i < samples.size(); ++i){
    index = (int)samples[i];

    for(int j = 0; j < x.ncol(); ++j){
      tmp = toUnitScale(x(index,j), lower, upper);
      rv = (R::runif(0,1)-0.5)*step;
      tmp = tmp + rv;
      if(tmp < 0) tmp = -tmp;
      if(tmp > 1) tmp = 2-tmp;
      x(index,j) = fromUnitScale(tmp, lower, upper);
    }
  }
  return(x);
} 
//[[Rcpp::export]]
NumericMatrix infection_proposal_function(NumericMatrix current, double lower, double upper, NumericMatrix step, int N, NumericMatrix iter, NumericVector indices){
  double move,tmp;
  NumericVector samples;
  NumericMatrix x(current.nrow(),current.ncol());
  int index = 0;
  int col = 0;
  //  samples = runif(N)*current.nrow();
  samples = randomSample(indices, N);
  for(int i = 0; i < samples.size(); ++i){
    index = (int)samples[i];
    col = (int)(R::runif(0,1)*3);
    tmp = toUnitScale(current(index,col),lower,upper);
    move = (R::runif(0,1)-0.5)*step(index,col);
    //move = R::rnorm(0.0,step(index,col)*step(index,col));
    tmp = tmp + move;
    if(tmp < 0) tmp = -tmp;
    if(tmp > 1) tmp = 2-tmp;
    tmp = fromUnitScale(tmp,lower,upper);
    current(index,col) = tmp;
    iter(index,col)++;
    x(index,col)++;
  }    
  return(x);
}

//[[Rcpp::export]]
double testings(){
  return(R::runif(-1,1));
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
double scaletuninglinear(double step, double popt, double pcur){
  double change = R::qnorm(popt/2,0,1,1,0)/R::qnorm(pcur/2,0,1,1,0);
  if(pcur >=1) pcur=0.99;
  if(pcur <=0) pcur=0.01;
  return(step*=change);
}

//[[Rcpp::export]]
int scale_infection_samples(int cur, double popt, double pcur, int N){
  double change = R::qnorm(popt/2,0,1,1,0)/R::qnorm(pcur/2,0,1,1,0);
  if(cur >= N) cur = N-1;
  if(cur <= 1) cur = 1;
  if(pcur >= 1) pcur=0.99;
  if(pcur<=0) pcur = 0.01;
  cur = cur*change;
  return(cur);
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
		double final_time,
		NumericMatrix actual_infection_times,
		NumericMatrix infection_steps
		){
  int infection_samples = 25; // Adaptive parameter dictating how many infection times are updated in one proposal step
  ofstream csv_write; // Output stream for csv results
  double TUNING_ERROR = 0.1; // Allowable tuning error for adaptive step sizes
  double newprobab, probab, difflike, single_weight,samp_weight; // Doubles for log-likelihoods and weights
  int no_recorded = 0; // Record number of recorded iterations
  int sampno = 1; // Record current sample number (note this may be different to no_recorded if thin != 1)
  int j = 0; // Record current parameter index
  string mcmc_chain_file = filename + "_chain.csv";
  string infection_times_file = filename + "_infectiontimes.csv";
  
  bool UPDATE_INFECTION_TIMES = false; // If true, updates infection times on step. If false, update other parameters

  // Numeric vectors for boost/waning parameters
  NumericVector startvalue(param_table.nrow()-1); // Starting values
  NumericVector step_sizes(param_table.nrow()-1); // Step sizes for proposal
  NumericVector current_params(param_table.nrow()-1); // Store current parameters
  NumericVector sampnos(save_block); // Stores sample numbers
  NumericVector lnlikes(save_block); // Stores log likelihoods
  NumericVector proposal(param_table.nrow()); // Proposed update to boost/wane parameters

  // Vectors for adaptive steps
  NumericVector tempaccepted(param_table.nrow()); // Store total number of accepted proposals for each parameter
  NumericVector tempiter(param_table.nrow()); // Store total number of proposals for each parameter
  NumericVector reset(param_table.nrow()); // Vector of zeroes for reset of above vectors
  NumericVector pcur(param_table.nrow()); // Vector of acceptance rates for each parameter

  NumericVector lower_bounds(param_table.nrow()); // Lower bounds for each parameter
  NumericVector upper_bounds(param_table.nrow()); // Upper bounds for each parameter
  NumericVector weights(param_table.nrow()); // Relative sampling weights for all parameters - this will become cumulative

  // Matrices to store the parameter set at each time point
  NumericMatrix empty_chain(save_block, (param_table.nrow()-1));
  NumericMatrix chain(save_block, (param_table.nrow()-1));

  // Matrix for proposed infection tmie update step
  NumericMatrix proposal_infection_times(infection_times.ncol(),infection_times.nrow());

  // Step sizes for infection time updates
  NumericMatrix ti_step_sizes(infection_steps.nrow(),infection_steps.ncol());;

  // Stores number of proposals and accepted proposals for ti parameters, as well as acceptance probability
  NumericMatrix ti_accepted(infection_steps.nrow(),infection_steps.ncol());
  NumericMatrix ti_iter(infection_steps.nrow(),infection_steps.ncol());
  NumericMatrix ti_pcur(infection_steps.nrow(),infection_steps.ncol()); // What if no proposals?
  NumericMatrix ti_reset(infection_steps.nrow(),infection_steps.ncol());
  NumericMatrix index_record(infection_steps.nrow(),infection_steps.ncol());

  NumericVector sample_indices(infection_steps.nrow());
  
  for(int i = 0; i < sample_indices.size();++i){
    sample_indices[i] = i;
  }


  // Really for debugging purposes - stores the differences between the current infection times and the actual infection times
  NumericVector time_differences(iterations + adaptive_period + burnin);
  
  // Take a copy of the given initial step sizes
  ti_step_sizes = clone(infection_steps);

  // Copy the given initial infection times to begin
  proposal_infection_times = clone(infection_times);
  
  // Get control parameters from parameter table
  for(int i = 0; i < (param_table.nrow()-1); ++i){
    startvalue[i] = param_table(i,1);
    step_sizes[i] = param_table(i,4);
  }
  lower_bounds = param_table(_,2);
  upper_bounds = param_table(_,3);


  // Convert given weights to scale 0-1 for sampling
  single_weight = 1/(accumulate(param_table(_,5).begin(), param_table(_,5).end(), 0.0));

  if(VERBOSE) cout << "Weights: ";
  weights[0] = param_table(0,5) * single_weight;
  if(VERBOSE)  cout << weights[0] << " ";
  for(int f = 1; f < (param_table.nrow());++f){
    weights[f] = param_table(f,5)*single_weight + weights[f-1];
    if(VERBOSE) cout << weights[f] << " ";
  }
  cout << endl;

  // Open stream to csv file
  csv_write.open(mcmc_chain_file.c_str());
  
  // Get initial parameter values and likelihood
  current_params = clone(startvalue);
  probab = posterior(current_params, data, proposal_infection_times, t0_titres, final_time);
  
  // Write first line to csv file
  csv_write << sampno << ",";
  if(VERBOSE) cout << "Starting parameters: ";
  for(int i = 0; i < (startvalue.size()); ++i){
    csv_write << startvalue[i] << ",";
    if(VERBOSE) cout << startvalue[i] << " ";
  }
  csv_write << probab << endl;
  if(VERBOSE) cout << endl;
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
    // Choose whether to update a parameter or infection times
    if(j == (weights.size()-1)) {
      index_record = clone(ti_reset);
      proposal_infection_times = clone(infection_times); // Copy current infection times

      // Propose update to infection times
      index_record = infection_proposal_function(proposal_infection_times, lower_bounds(j), upper_bounds(j), ti_step_sizes,infection_samples,ti_iter,sample_indices);

      //New likelihood
      newprobab = posterior(current_params, data, proposal_infection_times, t0_titres, final_time);

      //Difference
      difflike = newprobab - probab;
      
      //Metropolis acceptance
      if(R::runif(0,1) < exp(difflike) || difflike > 0){

	//	cout << difflike << endl;
	double tmp_sum = 0;
	for(int iii =0;iii<actual_infection_times.nrow();++iii){
	  for(int jjj=0;jjj<actual_infection_times.ncol();++jjj){
	    tmp_sum += abs(actual_infection_times(iii,jjj) - proposal_infection_times(iii,jjj));
	  }
	}
	time_differences[ii] = tmp_sum;
	
        infection_times = clone(proposal_infection_times);
	probab = newprobab;
	
	// Update acceptances
	for(int iiii = 0; iiii < ti_accepted.nrow();++iiii){
	  for(int jjjj=0;jjjj<ti_accepted.ncol();++jjjj){
	    ti_accepted(iiii,jjjj) += index_record(iiii,jjjj);
	  }
	}
	tempaccepted(j) += 1;
      }
      tempiter(j) += 1;
      //cout << "Iterations: " << tempiter(j) << endl;
      //cout << "Accepted: " << tempaccepted(j) << endl;
    }
    else{
      // Propose new value for current parameter
      proposal = clone(current_params);
      proposal(j) = proposalfunction(proposal(j),lower_bounds(j),upper_bounds(j),step_sizes(j));
      //proposal(j) = proposal_function(proposal(j),step_sizes(j));
      
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
      
    }
    
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
    if(opt_freq != 0 && ii < adaptive_period && (ii+1)%opt_freq==0){

          // Update step sizes for normal parameters
      for(int x = 0; x < current_params.size();++x){
	// get current acceptance rate
	if(tempiter(x) > 0){
	  pcur[x] = tempaccepted(x)/tempiter(x);
	  // If outside acceptable range:
	  if(pcur(x) < popt - (TUNING_ERROR*popt) || pcur(x) > popt + (TUNING_ERROR*popt)){
	    step_sizes(x) = scaletuning(step_sizes(x), popt, pcur(x));
	  }
	}
      }
      
      // 
      pcur[pcur.size()-1] = tempaccepted[tempaccepted.size()-1]/tempiter[tempiter.size()-1];
      /* if(pcur(pcur.size()-1) < popt - (TUNING_ERROR*popt) || pcur(pcur.size()-1) > popt + (TUNING_ERROR*popt)){
	 infection_samples = scale_infection_samples(infection_samples, popt, pcur(pcur.size()-1),infection_times.nrow());
	 }*/
      
      // Update step sizes of infection times
      for(int y = 0; y < ti_accepted.nrow();++y){
	for(int x=0; x < ti_accepted.ncol();++x){ 
	  if(ti_iter(y,x) > 0){
	    ti_pcur(y,x) = ti_accepted(y,x)/ti_iter(y,x);
	    if(ti_pcur(y,x) < popt - (TUNING_ERROR*popt) || ti_pcur(y,x) > popt + (TUNING_ERROR*popt)){
	      ti_step_sizes(y,x) = scaletuning(ti_step_sizes(y,x),popt,ti_pcur(y,x));
	    }
	    //	    cout << ti_pcur(y,x) << " ";;
	  }
	}
      }
      ti_accepted = clone(ti_reset);
      ti_iter = clone(ti_reset);
      tempaccepted = clone(reset);
      tempiter = clone(reset);
    }
  }
  csv_write.close();

  // Write infection times to csv
  csv_write.open(infection_times_file.c_str());
  for(int i = 0; i < infection_times.nrow(); ++i){
    for(int j = 0; j < infection_times.ncol(); ++j){
      csv_write << infection_times(i,j) << ",";
    }
    csv_write << endl;
  }
  csv_write.close();

  string time_differences_file = "time_differences.csv";
  csv_write.open(time_differences_file.c_str());
  for(int i = 0; i < time_differences.size();++i){
    csv_write << time_differences[i] << endl;
  }
  csv_write.close();

  return(mcmc_chain_file);
}
