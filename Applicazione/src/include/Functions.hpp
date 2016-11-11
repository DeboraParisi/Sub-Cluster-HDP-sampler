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
 *  In questo file sono raccolte tutte quelle funzioni di supporto all'algoritmo che non sono specifiche del modello scelto.
 *  \date Febbraio 2016
 */

 /*! \brief Criterio di confronto tra due elementi che sono coppie (unsigned int, double), da utilizzare nel sort di un vettore;
  *   il confronto è sul secondo valore nella coppia.
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
/*! \biref Permette di calcolare i numeri di stirling
 *	\param x1 primo elemento
 * 	\param x2 secondo elemento
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


/*! \brief Calcola i numeri di Stirling di prima specie s(n,m), per \f$ n = 0, \dots ,N \f$ e li memorizza nel vettore in ingresso.
 *  Per definizione: s(0,0) = s(1,1) = 1, s(n,0) = 0 per n > 0, s(n,m) = 0 per m > n, s(n,m) = s(n-1,m-1) + (n-1)*s(n-1,m).
 *  Il calcolo è in scala logaritmica per ottenere maggior precisione.
 *  \param N - valore massimo di n
 *  \param logstirling - vettore in cui memorizzare i numeri di Stirling in scala logaritmica
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

/*! \brief Kahan_algorithm
 *  Fa le somme di vettori e riduce l’errore numerico, dovuto all’arrotondamento che fa la macchina
 *  \param numbers - elementi da sommare tra di loro
 *  \return somma
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

/*! \brief Antoniak
 *  si occupa del campionamento dei tavoli in ogni gruppo, equazione 3.3 Relazione_Parisi_Perego.pdf
 *  \param Alpha - parametro di concentrazione del processo di dirchlet che governa i cluster nel gruppo
 *  \param Beta - peso globale del cluster k di cui vogliamo campionare i tavoli
 *  \param njk - numero di elementi del gruppo j nel cluster k
 *  \param LogStirling - Vettore dei numeri di Stirling in scala logaritmica
 *  \param Gen - Generatore di numeri casuali
 *  \return numero di tavoli che nel gruppo j servono il piatto k
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

/*! \brief Campionare gli \f$ \tilde{m_{jb}} \f$ , ovvero il numero di tavoli che in un documento j hanno il cluster b
 *  campiona i tavoli temporanei dal'equazione 3.24 Relazione_Parisi_Perego.pdf
 *  \param alpha - parametro di concentrazione del processo di dirchlet che governa i cluster nel gruppo
 *  \param K - numero corrente dei cluster
 *  \param njk - numero di elementi del gruppo j nel cluster k
 *  \param StirlingNumber - Vettore dei numeri di Stirling in scala logaritmica
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
