#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>
#include "omprng.hpp"
#include <math.h>   //funzione gamma
#include <random> // per campionare dalla uniforme
#include <algorithm> //per sort
#include <iostream>
#include <cmath>
#include <utility>
#include <iomanip>
#include <omp.h>
#include <limits>

using std::vector;

/*! \file Functions.hpp
 *  This file gathers useful functions for the algorithm, that are not model specific.
 *  \authors{Debora Parisi and Stefania Perego}
 *  \date February 2016
 */

/*! \brief Sorting operator for <unsigned int,double> pairs.
*    The order is based on the second element in the pair.
*    \authors{Debora Parisi and Stefania Perego}
*	 \date February 2016
*/
struct greater_for_pair{

       bool operator() (const std::pair<unsigned int,double>& x,const std::pair<unsigned int,double>& y)const {return x.second > y.second;}
};


/*function KahanSum(input)
    var sum = 0.0
    var c = 0.0                  // A running compensation for lost low-order bits.
    for i = 1 to input.length do
        var y = input[i] - c     // So far, so good: c is zero.
        var t = sum + y          // Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y        // (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
        sum = t                  // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
                                 // Next time around, the lost low part will be added to y in a fresh attempt.
    return sum
*/
/*! \brief Allows to compute Stirling numbers
 *	\param x1 first element
 * 	\param x2 second element
 */
long double logsumexp(long double x1,long double x2){

    if(x1>x2)
        return x1 + std::log(1 + std::exp(x2-x1));
    else if (x2>x1)
        return x2 + std::log(1 + std::exp(x1-x2));
    else if (x1 == std::numeric_limits<long double>::infinity())
        return std::numeric_limits<long double>::infinity();
    else if (x1 == -std::numeric_limits<long double>::infinity())
        return -std::numeric_limits<long double>::infinity();
    else
        return x1 + std::log(2);

}


/*! \brief Computes unsigned Stirling numbers of the first kind |s(n,m)|, for \f$ n = 0, \dots ,N \f$, and stores them in the input vector.
 *  By definition: s(0,0) = s(1,1) = 1, s(n,0) = 0 for n > 0, s(n,m) = 0 for m > n, s(n,m) = s(n-1,m-1) + (n-1)*s(n-1,m).
 *  To achieve greater precision, computes the numbers' logarithm.
 *  \param N - max value for n
 *  \param logstirling - vector that stores the logarithm of Stirling numbers
 *  \authors{Debora Parisi and Stefania Perego}
 *  \date February 2016
 */
void ComputeLogStirlingNumbers (unsigned int N, vector<long double> & logstirling){
	//riga 0
	logstirling[0]=0;

         if(N==0)
            return;

	 //riga 1
	 logstirling[1]= - std::numeric_limits<long double>::infinity();
	 logstirling[2]=0;

         if(N==1)
            return;

	 for (size_t k=1; k<N; ++k){
	     //ciclo for interno, fisso la riga e scorro le colonne della matrice
	     //fisso la riga n, devo scorrere stirling dalla posizione K a K+n
	     unsigned int n= k*(k+1)/2;
	     unsigned int n_plus_1 = (k+1)*(k+2)/2;

	     //stirling[n_plus_1]=0;
		 logstirling[n_plus_1]= -std::numeric_limits<long double>::infinity();

	     //stirling[n_plus_1 + (k+1)]=1;
		 logstirling[n_plus_1 + (k+1)]=0;

	     for (size_t i=1; i<(k+1); ++i)
			  logstirling[n_plus_1+i]= logsumexp( logstirling[n+i-1], + std::log(k) + logstirling[n+i]);
    }
}

/*! \brief Kahan_algorithm: computes the sum of the numbers contained in the input vector and reduces the round-off error due to the machine.  
 *  \param numbers - numbers to be summed 
 *  \return sum
 *  \authors{Debora Parisi and Stefania Perego}
 *  \date February 2016
 */
template<typename T> T Kahan_algorithm (vector<T>& numbers){

    vector<T> Temp_numbers = numbers;
    T sum = 0.0;
    T c = 0.0; // compensazione per la perdita degli ordini significativi più bassi
    //prima di iniziare ordinamo gli elementi del vettore in ordine crescente
    T y , t;
    std::sort (Temp_numbers.begin(), Temp_numbers.end());

    for(typename vector<T>::size_type i = 0; i < Temp_numbers.size(); i++){

        y= Temp_numbers[i] - c;
        t = sum + y;
        c = (t - sum) -y;
        sum = t;
    }

    return sum;
}

/*! \brief Antoniak: samples tables in each group. 
 *  \param Alpha - concentration parameter of the Dirichlet process governing a group
 *  \param Beta - global weight for cluster k
 *  \param njk - number of element of group j in cluster k
 *  \param LogStirling - vector of the logarithm of the Stirling numbers
 *  \param Gen - parallel random number generator
 *  \return number of tables in group j serving dish k
 *  \authors{Debora Parisi and Stefania Perego}
 *  \date February 2016
 */
unsigned int Antoniak (double alpha, double beta, unsigned int njk, const std::vector<long double>& LogStirling,omprng& Gen){

    //dato njk, mi servono gli stirlign s(njk,m) con m da 0 a njk e poi devo calcolare i pesi corrispondenti

    unsigned int m=0; //valore di ritorno
	double log_gamma_num= lgamma(alpha*beta);
	double log_gamma_den= lgamma(alpha*beta+njk);
	double alpha_beta= alpha*beta;
    unsigned int k = njk* (njk+1)/2;
    auto it= LogStirling.begin()+ k;
	vector<std::pair<unsigned int,double>> weights;
	weights.reserve(njk+1);
	double log_temp_weight=0.0;
	vector<double> TempForSum;

	if((njk==0)){              // m da 0 a njk, se njk !=0 m=0 ha prob nulla, ma se njk=0 allora m=0 con prob=1
		weights.push_back(std::make_pair(0,1.0));
        return 0;
    }
	else{
        weights.push_back(std::make_pair(0,0.0));
		for(size_t i=1; i< njk+1;i++){
            log_temp_weight= log_gamma_num-log_gamma_den+(*(it+i))+i*std::log(alpha_beta);
			weights.push_back(std::make_pair(i,std::exp(log_temp_weight)));
			}
		}

	double sum1=0.0;

	for(auto it=weights.cbegin(); it!=weights.cend(); it++)
	    TempForSum.push_back((*(it)).second);

	sum1 = Kahan_algorithm(TempForSum);

    for(size_t i=0; i<njk+1; i++)
        weights[i].second =(weights[i].second)/sum1;

	TempForSum.clear();

	for(auto it=weights.cbegin(); it!=weights.cend(); it++)
		TempForSum.push_back((*(it)).second);


	sum1 = Kahan_algorithm(TempForSum);

	if( !((std::fabs(1.0 - sum1)/std::fabs(1.0)) < 1.16e-09)){
        std::cerr<<"Error in Antoniak"<<std::endl;
        std::cerr<<"weights do no sum up to 1 "<<std::endl;
        exit(1);
	}

	//ordino i pesi in senso decrescente
	sort( weights.begin(), weights.end(), greater_for_pair() );

	double u = Gen.runif();

	// campiono m

	double weight_sum = 0.0;
	bool check = false;

    for(size_t i=0; i<njk+1;i++){
        weight_sum += weights[i].second;
        if( u <= weight_sum ){
			m = weights[i].first;
            check=true;
            break;
        }
    }

    if(!check){
		std::cerr<<"Error in Antoniak, no sample"<<std::endl;
		exit(1);
    }

    return m;
}

/*! \brief Computes \f$ \tilde{m_{jb}} \f$ , which approximates the tables in document j serving dish b 
 *  \param alpha - concentration parameter of the Dirichlet process governing a group
 *  \param K - current number of clusters
 *  \param njk - number of data of groupp j in cluster k
 *  \param StirlingNumber - vector of the logarithm of the Stirling numbers
 *  \authors{Debora Parisi and Stefania Perego}
 *  \date February 2016
*/
unsigned int FindBestNumTable(double alpha, unsigned int K,unsigned int njk, vector<long double> &StirlingNumber){

    std::pair <unsigned int, double> m (0, - std::numeric_limits<long double>::infinity());
    long double ArgMax;
    double gamma_Alpha_K=lgamma(alpha/K);
    long double gamma_Alpha_K_njk=lgamma((alpha/K) + njk);
    // mi sto posizionando nella riga di interesse, utilizzo l'iteratore
	unsigned int k = njk* (njk+1)/2;
	auto it=StirlingNumber.begin() + k;

    for(size_t i=0; i<=njk; i++){
        ArgMax=gamma_Alpha_K - gamma_Alpha_K_njk + (*(it+i)) + i*(std::log(alpha)-std::log(K));
        if(ArgMax > m.second){
            m.first=i;
            m.second=ArgMax;
        }
    }
    return m.first;
}

#endif
