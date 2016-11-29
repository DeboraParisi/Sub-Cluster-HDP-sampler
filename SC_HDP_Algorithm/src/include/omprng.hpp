#ifndef OMPRNG_HPP
#define OMPRNG_HPP

#include <omp.h>
#include <iostream>
#include "rngstream.hpp"
#include "sys/time.h"
#include <cmath>
#include <vector>

/*! \file omprng.hpp
 *
 * \brief Parallel random number generator for OpenMP
 */

/*! \brief Unnamed namespace where some constants are saved
	*/
namespace
{
	const double ENepero = std::exp(1.0);

	const double pi = std::acos(-1.0);
}

using namespace std;

/*! \brief Omprng library for sampling random numbers in OpenMp 
 *
 * The library was developed by
 * Matthew Bognar
 * Department of Statistics and Actuarial Science
 * University of Iowa
 * http://www.stat.uiowa.edu/~mbognar/omprng
 * matthew-bognar@uiowa.edu
 */
class omprng
{
	private:

		/*! \brief Number of processors available
		 */
		int nprocs;

		/*! \brief RngStream object. Implemented by Matthew Bognar
		 */
		RngStream *myRng;

	public:

		/*! \brief DEfault constructor
		 */
		omprng ();

		/*! \brief Destructor
		 */
		~omprng(){ 	delete[] myRng; }

		/*! \brief Set the seed
		 *  \param Seed - seed
			*/
		void fixedSeed(long unsigned int);

		/*! \brief Generate a random seed
			*/
		void randomSeed ();

		/*! \brief Set the number of threads
		 *  \param NumThread - number of threads
			*/
		void setNumThreads (int);

		/*! \brief Samples from the Uniform distribution between 0 and 1
		 *  \return a sample from the Uniform distribution between 0 and 1
			*/
		double runif ();

		/*! \brief Samples from the Uniform distribution between two fixed values
		 *  \param a - lower bound
		 *  \param b - upper bound
		 *  \return a sample from the Uniform distribution between two fixed values
         */
		double runif (double,double);

		/*! \brief Samples from the Gaussian distribution
		 *  \param mu - mean
		 *  \param sigma - standard deviation
		 *  \return a sample from the Gaussian distribution with given mean and variance
        */
		double rnorm (double,double);

		/** \brief Samples from the Exponential distribution
		 *
		 * X~exp(theta)
		 *
		 * f(x) = 1 /theta * exp(-x/theta)
		 *
		 * x > 0, theta > 0
		 *
		 * E(X) = theta, Var(X) = theta^2
		 * \param theta - scale parameter
		 * \return a sample from the Exponential distribution with given scale parameter
        */
		double rexp (double);

		/*! \brief Samples from the Gamma distribution
		 *
		 * X~gamma(alpha,beta)
		 *
		 * f(x) = 1/(gamma(alpha)*beta^alpha) * x^(alpha-1) * exp(-x/beta)
		 *
		 * x > 0, alpha > 0, beta > 0
		 *
		 * E(X) = alpha*beta, Var(X) = alpha*beta^2
		 * \param alpha - shape parameter
		 * \param beta - scale parameter
		 * \return a sample from the Gamma distribution with given shape and scale parameters
			*/
		double rgamma (double,double);

		/*! \brief Samples from the Chi-squared distribution
		 *
		 * X~chisq(df)
		 *
		 * f(x) = 1/(gamma(df/2)*2^(df/2)) * x^(df/2-1) * exp(-x/2)
		 *
		 * x > 0, df = 1,2,3,...
		 *
		 * E(X) = df, Var(X) = 2*df
		 * \param df - degrees of freedom  (dof)
		 * \return a sample from the Chi-squared distribution with given dof
			*/
		double rchisq (unsigned int);

		/*! \brief Samples from the Beta distribution
		 *
		 * X~beta(alpha,beta)
		 *
		 * f(x) = gamma(alpha+beta)/(gamma(alpha)*gamma(beta)) * x^(alpha-1) * (1-x)^(beta-1)
		 *
		 * 0 < x < 1, alpha > 0, beta > 0
		 *
		 * E(X) = alpha/(alpha+beta), Var(X) = alpha*beta / ((alpha+beta+1)*(alpha+beta)^2)
		 * \param alpha - first shape parameter
		 * \param beta - second shape parameter
		 * \return x - a sample from the Beta distribution with given shape parameters
			*/
		double rbeta (double,double);

		/*! \brief  Samples a discrete random variable with support 0:(K-1)
		 *
		 * K = size(inputvector)
		 * \param Logp - logarithm of weights associated to atoms
		 * \return a sample of a discrete random variable with given support ad weights
			*/
		unsigned int rdiscrete(std::vector<double> &);

		/*! \brief  Samples a discrete uniform random variable with support 0:(N-1) 
		 * \param N - number of classes
		 * \return x - a sample of a discrete uniform random variable with support 0:(N-1) 
         */
		unsigned int runifdiscrete(unsigned int);

		/*! \brief Samples from the Bernoulii distribution
		 *  X~Bernoulli(p)
		 *  P(X=1) = p
		 *  P(X=0) = 1-p
		 *  X={0,1}
		 *  \param p - success probability
		 *  \return a sample from the Bernoulii distribution with given probability of success
		 *  \authors{Debora Parisi and Stefania Perego}
		 *  \date February 2016
		 */
		unsigned int rbernoulli(double p);

        /*! \brief Samples from the Binomial distribution
		 *  X~Bin(n,p)
		 *
		 *  P(X=k) = n!/(k!(n-k)!) p^k (1-p)^(n-k)
		 *  k={0,1,...,n}
		 *  \param n - number of trials
		 *  \param p - success probability
		 *  \return a sample from the Binomial distribution with given number of trials and probability of success
		 *  \authors{Debora Parisi and Stefania Perego}
		 *  \date February 2016
		 */
		unsigned int rbinomial(unsigned int n,double p);

        /*! \brief Samples from the d-dimensional Dirichlet distribution
		 *  (X_1,...,X_d)~Dir(a_1,...,a_d)
		 *  a_i >0 i={1,..,d}
		 *  sum_{X_i} = 1
		 *  \param params - vector of parameters
         *  \param dir_sampled - a sample from the d-dimensional Dirichlet distribution with given parameters
		 *  \authors{Debora Parisi and Stefania Perego}
		 *  \date February 2016
		 */
		void rdirichlet(const vector<double>& params, vector<double>& dir_sampled);
};


#endif
