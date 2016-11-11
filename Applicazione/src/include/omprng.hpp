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
 * \brief Generatore di numeri casuali in parallelo per OpenMP
 */

/*! \brief Unnamed namespace in cui sono memorizzati delle costanti
	*/
namespace
{
	const double ENepero = std::exp(1.0);

	const double pi = std::acos(-1.0);
}

using namespace std;

/*! \brief Libreria Omprng per la generazione dei numeri casuali in OpenMp
 *
 * Si deve lo sviluppo a
 * Matthew Bognar
 * Department of Statistics and Actuarial Science
 * University of Iowa
 * http://www.stat.uiowa.edu/~mbognar/omprng
 * matthew-bognar@uiowa.edu
 * \date Luglio 2014
 */
class omprng
{
	private:

		/*! \brief Numero di processori a disposizione
		 */
		int nprocs;

		/*! \brief Oggetto RngStream. Implementazione dovuta a Matthew Bognar
		 */
		RngStream *myRng;

	public:

		/*! \brief Costruttore di default
		 */
		omprng ();

		/*! \brief Distruttore
		 */
		~omprng(){ 	delete[] myRng; }

		/*! \brief Fissa seed
		 *  \param Seed - seme casuale
			*/
		void fixedSeed(long unsigned int);

		/*! \brief Genera casualmente un seed
			*/
		void randomSeed ();

		/*! \brief Imposta il numero di thread
		 *  \param NumThread - numero di thread
			*/
		void setNumThreads (int);

		/*! \brief Generazione da Uniforme tra 0 ed 1
		 *  \return realizzazione di una uniforme tra 0 ed 1
			*/
		double runif ();

		/*! \brief Generazione da Uniforme tra due valori fissati
		 *  \param a - estremo inferiore
		 *  \param b - estremo superiore
		 *  \return realizzazione  di una uniforme tra a e b
         */
		double runif (double,double);

		/*! \brief Generazione da Gaussiana
		 *  \param mu - media
		 *  \param sigma - deviazione standard
		 *  \return realizzazione di una gaussiana con media mu e varianza sigma^2
        */
		double rnorm (double,double);

		/** \brief Generazione da Esponenziale.
		 *
		 * X~exp(theta)
		 *
		 * f(x) = 1 /theta * exp(-x/theta)
		 *
		 * x > 0, theta > 0
		 *
		 * E(X) = theta, Var(X) = theta^2
		 * \param theta - parametro di scala
		 * \return realizzazione di una esponenziale di parametro di scala theta
        */
		double rexp (double);

		/*! \brief Generazione da Gamma
		 *
		 * X~gamma(alpha,beta)
		 *
		 * f(x) = 1/(gamma(alpha)*beta^alpha) * x^(alpha-1) * exp(-x/beta)
		 *
		 * x > 0, alpha > 0, beta > 0
		 *
		 * E(X) = alpha*beta, Var(X) = alpha*beta^2
		 * \param alpha - parametro di forma
		 * \param beta - parametro di scala
		 * \return realizzazione di una gamma di parametro di forma alpha e di parametro di scala beta
			*/
		double rgamma (double,double);

		/*! \brief Generazione da ChiQuadro.
		 *
		 * X~chisq(df)
		 *
		 * f(x) = 1/(gamma(df/2)*2^(df/2)) * x^(df/2-1) * exp(-x/2)
		 *
		 * x > 0, df = 1,2,3,...
		 *
		 * E(X) = df, Var(X) = 2*df
		 * \param df - gradi di libertà
		 * \return realizzazione di una ChiQuadro con df gradi di libertà
			*/
		double rchisq (unsigned int);

		/*! \brief Generazione da Beta.
		 *
		 * X~beta(alpha,beta)
		 *
		 * f(x) = gamma(alpha+beta)/(gamma(alpha)*gamma(beta)) * x^(alpha-1) * (1-x)^(beta-1)
		 *
		 * 0 < x < 1, alpha > 0, beta > 0
		 *
		 * E(X) = alpha/(alpha+beta), Var(X) = alpha*beta / ((alpha+beta+1)*(alpha+beta)^2)
		 * \param alpha - primo parametro di forma
		 * \param beta - secondo parametro di forma
		 * \return x - realizzazione di una Beta con df parametri di forma alpha e beta
			*/
		double rbeta (double,double);

		/*! \brief Campionamento da una variabile aleatoria discreta con supporto 0:(K-1)
		 *
		 * K = size(inputvector)
		 * \param Logp - logaritmi probabilità/pesi
		 * \return realizzazione di una variabile aleatoria discreta con supporto 0:(K-1) e probabilità rispettive exp{Logp}
			*/
		unsigned int rdiscrete(std::vector<double> &);

		/*! \brief Campionamento da una variabile aleatoria discreta uniforme con supporto 0:(N-1)
		 * \param N - numero di classi per il campionamento discreto
		 * \return x - realizzazione di una variabile aleatoria discreta uniforme supporto 0:(N-1)
         */
		unsigned int runifdiscrete(unsigned int);

		/*! \brief Generazione da una Bernoulli
		 *  X~Bernoulli(p)
		 *  P(X=1) = p
		 *  P(X=0) = 1-p
		 *  X={0,1}
		 *  \param p - probabilità del successo
		 *  \return realizzazione di una Bernoulli con supporto {0,1}
		 *  \date Febbraio 2016
		 */
		unsigned int rbernoulli(double p);

        /*! \brief Generazione da una Binomiale
		 *  X~Bin(n,p)
		 *
		 *  P(X=k) = n!/(k!(n-k)!) p^k (1-p)^(n-k)
		 *  k={0,1,...,n}
		 *  \param n - numero di prove
		 *  \param p - probabilità di successo
		 *  \return realizzazione di una Binomiale con supporto {0,1,...,n}
		 *  \date Febbraio 2016
		 */
		unsigned int rbinomial(unsigned int n,double p);

        /*! \brief Generazione da una Dirichlet di dimensine d
		 *  (X_1,...,X_d)~Dir(a_1,...,a_d)
		 *  a_i >0 i={1,..,d}
		 *  sum_{X_i} = 1
		 *  \param params -  vettore dei parametri della Dirichlet
         *  \param dir_sampled - realizzazione della dirichlet
		 *  \date Febbraio 2016
		 */
		void rdirichlet(const vector<double>& params, vector<double>& dir_sampled);
};


#endif
