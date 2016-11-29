#ifndef _OMPRNG_CPP_
#define _OMPRNG_CPP_

#include "omprng.hpp"
#include "rngstream.hpp"
#include <omp.h>
#include <iostream>
#include "sys/time.h"
#include <cmath>
#include <vector>


omprng::omprng ()
{
	randomSeed();   				// default is random seed
	nprocs = omp_get_num_procs();	// determine number of processors
	myRng = new RngStream[nprocs];
	myRng -> AdvanceState(3,3);		// advance state
}

void omprng::fixedSeed(long unsigned int myInt) {
	unsigned long mySeed[6] = {myInt, myInt, myInt, myInt, myInt, myInt};
	RngStream::SetPackageSeed(mySeed);
	nprocs = omp_get_num_procs();
	delete[] myRng; 
	myRng = new RngStream[nprocs];
	myRng -> AdvanceState(3,3);
}

void omprng::randomSeed () {
	timeval tim;
	gettimeofday(&tim, NULL);
	unsigned long int seed[6] = {static_cast<unsigned long int>(tim.tv_sec),
								static_cast<unsigned long int>(tim.tv_usec),
								static_cast<unsigned long int>(tim.tv_sec+tim.tv_usec), 
		                         static_cast<unsigned long int>(abs(tim.tv_sec-tim.tv_usec)),
		                         static_cast<unsigned long int>(abs(tim.tv_usec-tim.tv_sec)), 5};
	RngStream::SetPackageSeed(seed);
}


void omprng::setNumThreads (int nt) {
	omp_set_num_threads(nt);		// set number of threads in OpenMP
	delete[] myRng;                // I dealoccate the memory
	myRng = new RngStream[nt];		// initialize RngStream
}

double omprng::runif () 
{
	//generate random number between 0 and 1, i.e. from Unif(0,1)

	return(myRng[omp_get_thread_num()].RandU01());
}

double omprng::runif (double a, double b)
{
	//generate random number between a and b, i.e. from Unif(a,b)

	return(a + (b-a)*myRng[omp_get_thread_num()].RandU01());
}

double omprng::rexp (double theta) 
{
	// generate random number from X~exp(theta)
	// f(x) = 1 /theta * exp(-x/theta)
	// x > 0, theta > 0
	// E(X) = theta, Var(X) = theta^2

	double u, x;

	u = runif();
	x = -log(u)*theta;

	return(x);
}

double omprng::rgamma (double alpha, double beta) 
{
	// generate random number from X~gamma(alpha,beta)
	// f(x) = 1/(gamma(alpha)*beta^alpha) * x^(alpha-1) * exp(-x/beta)
	// x > 0, alpha > 0, beta > 0
	// E(X) = alpha*beta, Var(X) = alpha*beta^2

	double x=0.0;
	double delta, v0, v1, v2, v3, psi=1.0, nu=1.0e10;

	for(int i=0; i<floor(alpha); i++)
		x = x + rexp(1.0);

	if(alpha > floor(alpha)) { // alpha has fractional part
		delta = alpha - floor(alpha);
		while(nu > pow(psi,delta-1.0)*exp(-psi)) {
			v1 = runif(); v2 = runif(); v3 = runif();
			v0 = ENepero / (ENepero + delta);
			if(v1 <= v0) {
				psi = pow(v2, 1.0/delta);
				nu = v3 * pow(psi,delta-1.0);
			}
			else {
				psi = 1.0 - log(v2);
				nu = v3 * exp(-psi);		
			}
		}
		x = x + psi;
	}

	return(beta*x);
}

double omprng::rbeta (double alpha, double beta)
{
	//generate random number from X~beta(alpha,beta)
	// f(x) = gamma(alpha+beta)/(gamma(alpha)*gamma(beta)) * x^(alpha-1) * (1-x)^(beta-1)
	// 0 < x < 1, alpha > 0, beta > 0
	// E(X) = alpha/(alpha+beta), Var(X) = alpha*beta / ((alpha+beta+1)*(alpha+beta)^2)

	double x, y;
	
	x = rgamma(alpha, 1.0);	
	
	if ( x == 0)
	   return 0.0;
	else{   
	   y = rgamma(beta, 1.0);
	   return(x/(x+y));
	}   
}

// Campionamento da una variabile aleatoria discreta con supporto 0:(K-1) 
// dove K = size(weights)

unsigned int omprng::rdiscrete(std::vector<double>& weights)
{

	double u = this -> runif();
	
	double sum_weights = 0;
	
	// Normalizzo i pesi passati in ingresso
	for(std::vector<double>::const_iterator it = weights.cbegin(); it!=weights.cend(); it++)
		sum_weights+=*it;
	
	double sum = 0;
	unsigned int sample = 0;
	
	for(std::vector<double>::const_iterator it = weights.cbegin(); it!=weights.cend(); it++)
    {
		sum += *it/sum_weights;

		if (u <= sum)
		{
			return sample;
		}
		
		sample++;
	}
	// Nel caso uniform sia pari ad 1
	return sample-1;
}


// Campionamento da una variabile aleatoria discreta uniforme con supporto 0:(N-1)
unsigned int omprng::runifdiscrete(unsigned int N)
{
		
	// Guess dell'uniforme
	double uniform = this -> runif(0,1);
	
	double sum = 0;
	double step=1.0/static_cast<double>(N);
	for(unsigned int n=0; n<N; n++)
    {
		sum += step;
		if (uniform <= sum)
		{
			return n;
		}
	}
	return N-1;
}

//Campionamento da una binomiale di parametri n,p: n nr prove p probabilità successo nella singola prova

unsigned int omprng::rbinomial(unsigned int n, double p){
  
    unsigned int x = 0;
	
	for(unsigned int i=0; i<n; ++i)
	    x += this -> rbernoulli(p);

    return x;
}

//Campionamento dalla bernoulli

unsigned int omprng::rbernoulli(double p){
    
	double u = this -> runif();
	
	if(u<p)
	  return 1;
	else
	  return 0;

}

//Campionamento dalla dirichlet 

void omprng::rdirichlet(const vector<double>& params, vector<double>& dir_sampled){
     
	unsigned int N = params.size(); 
	vector<double> gamma_sampled;
    gamma_sampled.reserve(N);
    double gamma_0=0.0; //la somma
    double sum = 0.0; // check
    double temp_sample= 0.0;

    for( std::vector<double>::size_type i=0; i<N; i++){
		temp_sample = this->rgamma(params[i],1.0);
		gamma_sampled.push_back(temp_sample);
	    gamma_0 += temp_sample;
    }
	
	//gamma_0 = Kahan_algorithm(gamma_sampled);
		
	if( gamma_0 == 0.0){
      std::cerr<< " Nella dirichlet non ho campionato perchè i parametri sono proprio piccoli"<<std::endl;	
	  exit(1);
	}else{

		for(size_t i=0; i<N; i++){
		   dir_sampled[i]= gamma_sampled[i] / gamma_0;
		   sum += dir_sampled[i];  
		}   

		//sum = Kahan_algorithm(dir_sampled);

		if( !((std::fabs(1.0 - sum)/std::fabs(1.0)) < 9.99e-10)){
            std::cerr<<"Error in sampling_dirichlet"<<std::endl;
            std::cerr<<"the sum of the component just sampled, it isn't 1 "<<std::endl;
            exit(1);
		}
    } // fine else

}

#endif 
