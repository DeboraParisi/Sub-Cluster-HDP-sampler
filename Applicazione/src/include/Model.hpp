#ifndef _MODEL_HPP_
#define _MODEL_HPP_

#include "Cluster.hpp"
#include "Functions.hpp"
#include "omprng.hpp"

#include <vector>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>

using std::vector;
using std::unordered_map;

/*! \file Model.hpp
 *
 * \brief Classe che si occupa di gestire tutti i cluster, campionamenti che dipendo dal modello scelto e
 * funzioni che dipendono dal modello
 * \date Febbraio 2016
 */

/*! \brief Modello Generico di Model
 *
 *  Classe astratta dove tutti i metodi virtuali sono null.
 *  Le classi che ereditano da ModelGeneric servono per campionare i parametri latenti e gestire i relativi iperparametri, che di pendono
 *  dal modello scelto ed utilizzare.
 *  Gestisce funzioni variano in base la modello, come le il calolo delle versimiglianze, delle marginali e delle densità.
 *  Gestione di tutti i cluster, aggiunge e rimuove cluster.
 * \date Febbraio 2016
 */

template<typename Type, unsigned int DIM=1>
class ModelGeneric{

    public:

        /*! \brief Da informazione sul numero di clusters correnti
         *  \return il numero corrente di clusters
         */
        virtual unsigned int ViewK() const = 0;

        /*! \brief Visualizza gli ID dei clusters attulmente presenti
         *  \param Key - oggetto dove salvare gli Id dei clusters
         */
        virtual void ViewKey(vector<unsigned int>& ) const = 0;

        /*! \brief Fissa gl iperparametr dei parametri latenti
         *  \param _Lambda - oggetto ti tipo HYP, iperparametri in ingresso del parametro latente
         */
        virtual void SetHyperparameter(const typename Type::HYP&) = 0;

        /*! \brief Fissa gli iperparametri dei parametri latenti con i valori di default
         *  \param W - dimensione degli iperparametri
         */
        virtual void DefaultHyperparameter(size_t) = 0;

        /*! \brief Fissa i cluster iniziali, assegnando i pesi iniziali
         *  \param _K - Numero di cluster iniziali
         */
        virtual void SetInitialClusters(unsigned int) = 0;

        /*! \brief Calcola la verosmiglianza logaritmica marginale del cluster, una volta specificato l'id del clusters
         *  \param _K - id del cluster di cui si vuole calcolare la verosmiglianza logaritmica marginale
         */
        virtual double Marginalized_Loglikelihood(const unsigned int)  = 0;

        /*! \brief Calcola la loglikelihood in un punto X una volta precisato l'id del cluster
         *  \param X - punto di valutazione
         *  \param _K - id del cluster
         *  \return verosmiglianza logaritmica valutante nel punto X del cluster _K
         */
        virtual double Loglikelihood(const typename Type::Point , const unsigned int) = 0;

        /*! \brief Calcola la loglikelihood in un punto X una volta precisato l'id del sub-cluster sinistro
         *  \param X - punto di valutazione
         *  \param _K - id del sub-cluster sinistro
         *  \return verosmiglianza logaritmica valutante nel punto X del sub-cluster sinistro del cluster _K
         */
        virtual double LoglikelihoodLeft(const typename Type::Point, const unsigned int) = 0;

        /*! \brief Calcola la loglikelihood in un punto X una volta precisato l'id del sub-cluster destro
         *  \param X - punto di valutazione
         *  \param _K - id del sub-cluster destro
         *  \return verosmiglianza logaritmica valutante nel punto X del sub-cluster destro del cluster _K
         */
        virtual double LoglikelihoodRight(const typename Type::Point, const unsigned int) = 0;

        /*! \brief Calcola la densità del parametro latente che caratterizza il cluster specificato.
         *  \param _K - Id del cluster
         *  \return Desità logaritmica del cluster _K
         */
        virtual long double LogDensity(const unsigned int) = 0;

        /*! \brief Aggirnamento dei parametro latenti di tutti i clusters.
         *
         *  Equazione (3.6) della relazione Relazione_Parisi_Perego.pdf
         *  \param Gen - Generatore dei numeri casuali
         */
        virtual void UpdateThetaCluster(omprng&) = 0;

        /*! \brief Aggirnamento dei parametro latenti del sub-cluster sinistro e destro di tutti i clusters.
         *
         *  Equazione (3.10) della relazione Relazione_Parisi_Perego.pdf
         *  \param Gen - Generatore dei numeri casuali
         */
        virtual void UpdateThetaSubCluster(omprng&) = 0;

        /*! \brief Aggirnamento dei parametro latenti di un cluster, una volta precisato il suo id.
         *
         *  Equazione (3.6) della relazione Relazione_Parisi_Perego.pdf
         *  \param _K - id del cluster di cui si vuole aggiornare parametri lantenti
         *  \param Gen - Generatore dei numeri casuali
         */
        virtual void UpdateOneThetaCluster(const unsigned int,omprng&) = 0;

        /*! \brief Aggirnamento dei parametro latenti del sub-cluster sinistro e destro di un cluster, una volta precisato il sui id.
         *
         *  Equazione (3.10) della relazione Relazione_Parisi_Perego.pdf
         *  \param _K - id del cluster di cui si vuole aggiornare i parametri latenti nei sui sub-clusters
         *  \param Gen - Generatore dei numeri casuali
         */
        virtual void UpdateOneThetaSubCluster(const unsigned int,omprng&) = 0;

        /*! \brief Aggiunge un cluster vuoto ai clusters già presenti ed aggiorna il parametro K, che identifica il numero di cluster attuali.
         *         I parametri del cluster devono essere specificati in un secondo momento
         *  \param _K - Id del nuovo cluster che si aggiunge.
         */
        virtual void AddOneCluster(const unsigned int) = 0;

        /*! \brief Rimuove un cluster una volta specificato l'id, ed aggiorna il parametro K, che identifica il numero di cluster attuali
         *  \param _K - Id del cluster che si vuole rimuovere
         */
        virtual void RemoveOneCluster (const unsigned int) = 0;

        /*! \brief Rimuove più clusters una volta che sono specificati gli ID, ed aggiorna il parametro K, che identifica il numero di cluster attuali
         *  \param _K - Vettore degli Id dei clusters che si vogliono rimuovere
         */
        virtual void RemoveClusters (const vector<unsigned int>&) = 0;

        /*! \brief Mostra i pesi globali dei cluster attualemente presenti
         *  \param AllBeta - vettore che verrà riempito con i pesi globali dei clusters
         */
        virtual void ViewBeta(vector<double>&) = 0;

        /*! \brief Mostra i pesi globali del sub-cluster sinistro dei clusters attualemente presenti
         *  \param AllBetaLeft - vettore che verrà riempito con i pesi globali del sub-cluster sinistro dei clusters
         */
        virtual void ViewBetaLeft(vector<double>&) = 0;

        /*! \brief Mostra i pesi globali del sub-cluster destro dei clusters attualemente presenti
         *  \param AllBetaRight - vettore che verrà riempito con i pesi globali del sub-cluster destro dei clusters
         */
        virtual void ViewBetaRight(vector<double>&) = 0;

        /*! \brief Stampa su file i valori dei parametri latenti dei cluster attualmente presenti
         */
        virtual void PrintTheta(const std::string&) = 0;

        /*! Stampa a video le informazioni degli iperparametri dei parametri latenti dei cluster
         */
        virtual void PrintLambdaInfo() const  = 0;

       // virtual void ViewModel() = 0;

};

/*! \brief Modello Categorical
 *
 *  Verosimiglianza: Categorical.
 *  Prior sui parametri latenti: Dirichlet
 *
 * Questa classe è impiegata per campionare i parametri latenti da una Dirichlet, per il calcolo delle verosimiglianze e delle marginali
 * per un modello Categorica.
 * Gestisce anche l'aggiunta e la rimozione dei clusters, la stampa su file dei parametri latenti.
 *
 *  \date Febbrario 2016
 */
template <unsigned int DIM=1>
class CategoricalModel final: public ModelGeneric<TypeCategorical<DIM>,DIM>{

    public:

        //DEFIINIZIONE DEI TIPI PUBBLICI


        /*! \brief Parametro latente, vettore di pesi degli elementi distinti nel cluster
         */
        using THETA = TypeCategorical<1>::THETA;

        /*! \brief Singlo dato, dati ripetuti
         */
        using POINT = TypeCategorical<1>::Point;

        /*! \brief Statistiche per aggiornare l'iperparametro del parametro latente, numero di dati che sono cotentuti nel cluster e sub-cluster
         */
        using STAT  = TypeCategorical<1>::STAT;

        /*! \brief Iperparametro del parametro latente
         */
        using HYP   = TypeCategorical<1>::HYP;

    private:

        //DEFINIZIONE DEI TIPI PRIVATI

        /*! \brief Vettore degli iperparametri
         */
        HYP Lambda;

        /*! \brief numero corrente di cluster
         */
        unsigned int K;

        /*! \brief Insieme degli oggetti di tipo CategoricaCluster, individuati in base al sui Id
         */
        unordered_map<unsigned int, CategoricalCluster> Clusters;

        /*! Numero di threads
         */
        unsigned int OMP_NUM_THREADS;

    public:

        //COSTRUTTORI E DISTRUTTORE

        /*! \brief Costruttore di default
         */
        CategoricalModel();

		/*! \brief Copy costructor
		 */
		CategoricalModel(const CategoricalModel& mod)= default;

		/*! \brief Move constructor
		 */
		CategoricalModel(const CategoricalModel&& mod);

		/*! \brief Distruttore di default
		 */
        ~CategoricalModel()=default;

        //OPERATORI

        /*! \brief Access operator per gli l-value
         *  \param _K - a quale cluster si vuole accedere
         */
        CategoricalCluster& operator[] (unsigned int _K);

		/*! \brief Assignement operator
		 *  \param mod - oggetto di tipo CategoricalModel
		 */
        CategoricalModel& operator=(const CategoricalModel &mod);

        /*! \brief Move assignement operator
		 *  \param mod - oggetto di tipo CategoricalModel
		 */
        CategoricalModel& operator=(CategoricalModel &&mod);

        /*! \brief Fissa gli iperparametri dei parametri latenti
         *  \param _Lambda - Iperparametri in ingresso
         */
        void SetHyperparameter(const HYP& _Lambda);

        /*! \brief Fissa gli iperparametri dei parametri latenti con i valori di default
         *  \param W - dimensione degli iperparametri
         */
        void DefaultHyperparameter(size_t W);

        /*! \brief Fissa i cluster iniziali, assegnando i pesi globali iniziali
         *  \param _K - Numero di cluster iniziali
         */
		void SetInitialClusters(unsigned int _K);

        //METODI SPECIFICI DEL MODELLO:
        //VERISIMIGLIANZE MARGINALE DENSITÀ

		/*! \brief Calcola la verosmiglianza logaritmica marginale del cluster, una volta specificato l'id del clusters
         *  \param _K - id del cluster di cui si vuole calcolare la verosmiglianza logaritmica marginale
         */
		double Marginalized_Loglikelihood(const unsigned int _K);

        /*! \brief Calcola la loglikelihood in un punto X una volta precisato l'id del cluster
         *  \param X - punto di valutazione
         *  \param _K - id del cluster
         *  \return verosmiglianza logaritmica valutante nel punto X del cluster _K
         */
        double Loglikelihood(const POINT X,const unsigned int _K);

        /*! \brief Calcola la loglikelihood in un punto X una volta precisato l'id del sub-cluster sinistro
         *  \param X - punto di valutazione
         *  \param _K - id del sub-cluster sinistro
         *  \return verosmiglianza logaritmica valutante nel punto X del sub-cluster sinistro del cluster _K
         */
        double LoglikelihoodLeft(const POINT X,const unsigned int _K);

        /*! \brief Calcola la loglikelihood in un punto X una volta precisato l'id del sub-cluster destro
         *  \param X - punto di valutazione
         *  \param _K - id del sub-cluster destro
         *  \return verosmiglianza logaritmica valutante nel punto X del sub-cluster destro del cluster _K
         */
        double LoglikelihoodRight(const POINT X,const unsigned int _K);

        /*! \brief Calcola la densità del parametro latente che caratterizza il cluster specificato.
         *  \param _K - Id del cluster
         *  \return Densità logaritmica del cluster _K
         */
        long double LogDensity(const unsigned int _K);

        //METODI PER AGGIORNARE I PARAMETRI LANTENTI

        /*! \brief Aggirnamento dei parametro latenti di tutti i clusters.
         *
         *  Equazione (3.6) della relazione Relazione_Parisi_Perego.pdf
         *  \param Gen - Generatore dei numeri casuali
         */
        void UpdateThetaCluster(omprng& Gen);

        /*! \brief Aggirnamento dei parametro latenti del sub-cluster sinistro e destro di tutti i clusters.
         *
         *  Equazione (3.10) della relazione Relazione_Parisi_Perego.pdf
         *  \param Gen - Generatore dei numeri casuali
         */
        void UpdateThetaSubCluster(omprng& Gen);

        /*! \brief Aggirnamento dei parametro latenti di un cluster, una volta precisato il suo id.
         *
         *  Equazione (3.6) della relazione Relazione_Parisi_Perego.pdf
         *  \param _K - id del cluster di cui si vuole aggiornare parametri lantenti
         *  \param Gen - Generatore dei numeri casuali
         */
        void UpdateOneThetaCluster(const unsigned int _K,omprng& Gen);

        /*! \brief Aggirnamento dei parametro latenti del sub-cluster sinistro e destro di un cluster, una volta precisato il sui id.
         *
         *  Equazione (3.10) della relazione Relazione_Parisi_Perego.pdf
         *  \param _K - id del cluster di cui si vuole aggiornare i parametri latenti nei sui sub-clusters
         *  \param Gen - Generatore dei numeri casuali
         */
        void UpdateOneThetaSubCluster(const unsigned int _K,omprng& Gen);

        //METODI PER AGGIUNGERE E RIMUOVERE I CLUSTER

        /*! \brief Aggiunge un cluster vuoto ai clusters già presenti ed aggiorna il parametro K, che identifica il numero di cluster attuali.
         *         I parametri del cluster devono essere specificati in un secondo momento
         *  \param _K - Id del nuovo cluster che si aggiunge.
         */
        void AddOneCluster (const unsigned int _k);

        /*! \brief Rimuove un cluster una volta specificato l'id, ed aggiorna il parametro K, che identifica il numero di cluster attuali
         *  \param _K - Id del cluster che si vuole rimuovere
         */
        void RemoveOneCluster(const unsigned int _K);

         /*! \brief Rimuove più clusters una volta che sono specificati gli ID, ed aggiorna il parametro K, che identifica il numero di cluster attuali
         *   \param _K - Vettore degli Id dei clusters che si vogliono rimuovere
         */
        void RemoveClusters(const vector<unsigned int>& _K);

        //ALTRI METODI
        /*! \brief Da informazione sul numero di clusters correnti
         *  \return il numero corrente di clusters
         */
        unsigned int ViewK() const;

        /*! \brief Visualizza gli ID dei clusters attulmente presenti
         *  \param Key - oggetto dove salvare gli Id dei clusters
         */
        void ViewKey(vector<unsigned int>& Key) const;

		/*! \brief Mostra i pesi globali dei cluster attualemente presenti
         *  \param AllBeta - vettore che verrà riempito con i pesi globali dei clusters
         */
	   	void ViewBeta(vector<double>& _AllBeta);

	   	/*! \brief Mostra i pesi globali del sub-cluster sinistro dei clusters attualemente presenti
         *  \param AllBetaLeft - vettore che verrà riempito con i pesi globali del sub-cluster sinistro dei clusters
         */
		void ViewBetaLeft(vector<double>& _AllBetaLeft) ;

		/*! \brief Mostra i pesi globali del sub-cluster destro dei clusters attualemente presenti
         *  \param AllBetaRight - vettore che verrà riempito con i pesi globali del sub-cluster destro dei clusters
         */
		void ViewBetaRight(vector<double>& _AllBetaRight);

        /*! \brief Stampa su file i valori dei parametri latenti dei cluster attualmente presenti, i pesi degli elementi distinti in ogni clusters
         */
		void PrintTheta(const std::string& FileName) ;

		/*! Stampa a video le informazioni degli iperparametri dei parametri latenti dei cluster
         */
		void PrintLambdaInfo() const;

		//void ViewModel();

};

///////////////////////////////////////
// DEFINIZIONI DI CATEGORICALMODEL ////
///////////////////////////////////////

//Costruttore di default
template<unsigned int DIM> CategoricalModel<DIM>::CategoricalModel(){
    if(DIM!=1){
        std::cerr<<"Error in CategoricalModel"<<std::endl;
        std::cerr<<"In CategoricalModel la dimensione puo' essere solo pari ad 1"<<std::endl;
        exit(1);
    }

	K = 0;

	// Valore di default numero massimo di thread (seriale)
	OMP_NUM_THREADS=1;
	// Lettura numero massimo di thread
	char* EnvironmentVariable;
	EnvironmentVariable = getenv("OMP_NUM_THREADS");
	if (EnvironmentVariable!=NULL)
	{
		std::stringstream ss;
		ss << EnvironmentVariable;
		ss >> OMP_NUM_THREADS;
	}


}

// Fissa i cluster iniziali, assegnando i pesi globali iniziali
template<unsigned int DIM> void CategoricalModel<DIM>::SetInitialClusters(unsigned int _K){

    CategoricalCluster New;

	for(size_t k=0; k<_K; ++k){
	    Clusters.insert({k,New});
		Clusters[k].SetBeta(1.0/(_K+1));
    }

    K = _K;

}


//move constructor
template<unsigned int DIM> CategoricalModel<DIM>::CategoricalModel(const CategoricalModel&& mod): Lambda(mod.Lambda),K(mod.K),Clusters(mod.Clusters){}

// copy assignment
template<unsigned int DIM> CategoricalModel<DIM>& CategoricalModel<DIM>::operator=(const CategoricalModel<DIM>& mod){

    Lambda = mod.Lambda;
	K = mod.K;
	Clusters = mod.Clusters;

	return *this;
}

// move assignment operator
template<unsigned int DIM> CategoricalModel<DIM>& CategoricalModel<DIM>::operator=(CategoricalModel<DIM>&& mod) {

    Lambda = mod.Lambda;
	K = mod.K;
	Clusters = mod.Clusters;

	return *this;
}

//Fissa gli iperparametri dei parametri latente
template<unsigned int DIM> void CategoricalModel<DIM>::SetHyperparameter(const HYP& _Lambda){

    for(size_t i=0; i<_Lambda.size(); i++){
        if(_Lambda[i]<0){
            std::cerr<<"Error in CategoricalModel, all elements of Lambda must be positive"<<std::endl;
            exit(1);
        }
    }

    Lambda=_Lambda;
}

//Fissa gli iperparametri dei parametri lantenti con dei valori di default
template<unsigned int DIM> void CategoricalModel<DIM>::DefaultHyperparameter(size_t W){

    Lambda.assign(W,1.0);
}

//Operatore di accesso ai cluster e ai loro metodi dall'esterno
template<unsigned int DIM> CategoricalCluster& CategoricalModel<DIM>::operator[] (unsigned int _K){
    return Clusters[_K];
}

//marginale dei dati
template <unsigned int DIM> double CategoricalModel<DIM>::Marginalized_Loglikelihood(const unsigned int _K)
{

	double m_ll = 0.0;
	double lambda_sum = Kahan_algorithm(Lambda);
	STAT cw;
	Clusters[_K].ViewStatistics(cw);
	unsigned int c = Kahan_algorithm(cw);

	m_ll += lgamma(lambda_sum) - lgamma(lambda_sum + c);

	for(size_t w=0; w<Lambda.size(); ++w)
		m_ll = lgamma(Lambda[w] + cw[w]) - lgamma(Lambda[w]);

	return m_ll;
}

//Verosimiglianza logaritmica valutata nel dato X una volta fissato il cluster _K
template <unsigned int DIM> double CategoricalModel<DIM>::Loglikelihood(const POINT X, const unsigned int _K){
    double log_temp=0;
    double theta_temp;

    theta_temp=Clusters[_K].ViewThetaId(X - 1);

    log_temp=std::log(theta_temp);

    return log_temp;
}

//Verosimiglianza logaritmica del sub-cluster sinistro valutata nel dato X una volta fissato il cluster _K
template <unsigned int DIM> double CategoricalModel<DIM>::LoglikelihoodLeft(const POINT X, const unsigned int _K){
    double log_temp=0;
    double theta_temp;

    theta_temp=Clusters[_K].ViewThetaLeftId(X - 1);

    log_temp=std::log(theta_temp);

    return log_temp;
}

//Verosimiglianza logaritmica del sub-cluster destro valutata nel dato X una volta fissato il cluster _K
template <unsigned int DIM> double CategoricalModel<DIM>::LoglikelihoodRight(const POINT X, const unsigned int _K){
    double log_temp=0;
    double theta_temp;

    theta_temp=Clusters[_K].ViewThetaRightId(X - 1);

    log_temp=std::log(theta_temp);

    return log_temp;
}

//Densità logaritmica
template <unsigned int DIM> long double CategoricalModel<DIM>::LogDensity(const unsigned int _K){

    long double LogD = 0.0;
    HYP param;
    STAT c_k;
    THETA Data;


    Clusters[_K].ViewStatistics(c_k);
    Clusters[_K].ViewTheta(Data);

    for(size_t i=0; i<c_k.size(); i++)
        param.push_back(c_k[i] + Lambda[i]);

    double Sum_Alpha = Kahan_algorithm(param);

    LogD += lgamma(Sum_Alpha);

	for(vector<double>::size_type n=0; n<param.size(); ++n){
		LogD+= - lgamma(param[n]) + (param[n] - 1) * log(Data[n]);
	}

	return LogD;

}
//Metodi di aggiornamento dei parametri theta
template<unsigned int DIM> void CategoricalModel<DIM>::UpdateThetaCluster(omprng& Gen){


    unsigned long N=Lambda.size();
    vector<double> params;
    STAT c;
    THETA temp_theta;


    Gen.setNumThreads(OMP_NUM_THREADS);

        #pragma omp parallel for private(c,params,temp_theta) shared(N) schedule(static,1)
        for(size_t k=0; k<K; k++){

		        Clusters[k].ViewStatistics(c);
     		    temp_theta.resize(N);
     		    params.resize(N);

                if(Lambda.size()!= c.size()){
                        std::cerr<<"Errore in CategoricalModel: UpdateThetaCluster"<<std::endl;
                        std::cerr<<"Lambda e il vettore dei contatori hanno dimensione differente"<<std::endl;
                        exit(1);
                }

				for(size_t i=0; i<N; i++)
                        params[i] = c[i]+Lambda[i];


                Gen.rdirichlet (params,temp_theta);
                Clusters[k].SetTheta(temp_theta);

                c.clear();
                params.clear();
	            temp_theta.clear();

                }

}

//Aggiorna i paramentri lantenti dei sub-cluster di tutti i clusters
template<unsigned int DIM> void CategoricalModel<DIM>::UpdateThetaSubCluster(omprng& Gen){

    vector<double> paramsLeft;
    vector<double> paramsRight;
    STAT cLeft;
    STAT cRight;
    unsigned long N=Lambda.size();
    THETA temp_thetaRight;
    THETA temp_thetaLeft;

    Gen.setNumThreads(OMP_NUM_THREADS);


	#pragma omp parallel for private (paramsLeft,paramsRight,cLeft,cRight,temp_thetaRight,temp_thetaLeft) shared(N) schedule(static,1)
    	for(size_t k=0; k<K; k++){

        Clusters[k].ViewStatisticsLeft(cLeft);
        paramsLeft.resize(N);
        Clusters[k].ViewStatisticsRight(cRight);
        paramsRight.resize(N);

        if(Lambda.size()!= cLeft.size()){
           	std::cerr<<"Errore in CategoricalModel: UpdateThetaCluster"<<std::endl;
            std::cerr<<"Lambda e il vettore dei contatori hanno dimensione differente"<<std::endl;
            exit(1);
        }
        if(Lambda.size()!=cRight.size()){
           	std::cerr<<"Errore in CategoricalModel: UpdateThetaCluster"<<std::endl;
            std::cerr<<"Lambda e il vettore dei contatori hanno dimensione differente"<<std::endl;
            exit(1);
        }
        for(size_t i=0; i<N; i++){
            paramsLeft[i] = cLeft[i]+Lambda[i];
            paramsRight[i] = cRight[i]+Lambda[i];
        }

        temp_thetaLeft.resize(N);
		Gen.rdirichlet (paramsLeft,temp_thetaLeft);

		temp_thetaRight.resize(N);
        Gen.rdirichlet(paramsRight,temp_thetaRight);

		Clusters[k].SetThetaLeft(temp_thetaLeft);
        Clusters[k].SetThetaRight(temp_thetaRight);
		cLeft.clear();
        paramsLeft.clear();
        cRight.clear();
        paramsRight.clear();


	}


}

//Aggiorna  parametri latenti di un cluster specifico, _K
template<unsigned int DIM> void CategoricalModel<DIM>::UpdateOneThetaCluster(const unsigned int _K,omprng& Gen){

    vector<double> params;
    STAT c;
    unsigned long N=Lambda.size();
    THETA temp_theta(N,0);

    Clusters[_K].ViewStatistics(c);

    if(Lambda.size()!= c.size()){
        std::cerr<<"Errore in CategoricalModel: UpdateOneThetaCluster"<<std::endl;
        std::cerr<<"Lambda e il vettore dei contatori hanno dimensione differente"<<std::endl;
        exit(1);
    }

	for(size_t i=0; i<N; i++){
        params.push_back(c[i]+Lambda[i]);
	}

    Gen.rdirichlet(params,temp_theta);
    Clusters[_K].SetTheta(temp_theta);

}

//Aggiorna i parametri lantenti dei sub-cluster di un cluster specifico, _K
template <unsigned int DIM> void CategoricalModel<DIM>::UpdateOneThetaSubCluster(const unsigned int _K, omprng& Gen){

    vector<double> paramsLeft;
    vector<double> paramsRight;
    STAT cLeft;
    STAT cRight;
    unsigned long N=Lambda.size();
    THETA temp_thetaRight(N, 0.0);
    THETA temp_thetaLeft(N, 0.0);


    Clusters[_K].ViewStatisticsLeft(cLeft);
    Clusters[_K].ViewStatisticsRight(cRight);

    if(Lambda.size()!= cLeft.size() && Lambda.size()!=cRight.size()){
        std::cerr<<"Errore in CategoricalModel: UpdateThetaCluster"<<std::endl;
        std::cerr<<"Lambda e il vettore dei contatori hanno dimensione differente"<<std::endl;
        exit(1);
    }

	for(size_t i=0; i<N; i++){
        paramsLeft.push_back(cLeft[i]+Lambda[i]);
        paramsRight.push_back(cRight[i]+Lambda[i]);
	}

	Gen.rdirichlet(paramsLeft,temp_thetaLeft);
	Gen.rdirichlet(paramsRight,temp_thetaRight);

    Clusters[_K].SetThetaLeft(temp_thetaLeft);
    Clusters[_K].SetThetaRight(temp_thetaRight);

}

//Aggiunge un cluster ed aggiorna il numero attuale di cluster K
template<unsigned int DIM> void CategoricalModel<DIM>::AddOneCluster(const unsigned int _k){

    CategoricalCluster temp;
    if(!Clusters.insert({_k,temp}).second){
        std::cerr<<"Errore in CategoricalModel: AddOneCluster"<<std::endl;
        std::cerr<<"Non riesco a generare una chiave non esistente"<<std::endl;
        exit(1);
    }

    K++;
}

//Rimuove un cluster specifico, _K. Aggiorna il numero attuale di cluster K
template<unsigned int DIM> void CategoricalModel<DIM>::RemoveOneCluster(const unsigned int _K){

    unordered_map<unsigned int, CategoricalCluster> ClustersNew;

	//inserisco gli elementi prima di quello da togliere (stesse chiavi)
	for(size_t k=0; k < _K ; ++k)
	    ClustersNew.insert(std::make_pair(k,Clusters[k]));
	//inserisco gli elementi dopo quello da togliere (chiave scalate all'indietro di 1)
	for(size_t i= (_K+1); i < Clusters.size() ; ++i)
	    ClustersNew.insert(std::make_pair(i-1,Clusters[i]));

    Clusters.swap(ClustersNew);

    K--;
}

//Rimuove un insieme di cluster di cui abbiamo specificato l'id, aggiorna poi il numero attuale di cluster presenti.
template<unsigned int DIM> void CategoricalModel<DIM>::RemoveClusters(const vector<unsigned int>& _K){

    if(_K.empty())
        return;

    unordered_map<unsigned int, CategoricalCluster> ClustersNew;
	bool IsPresent = false;
    unsigned int NewKey = 0;

   	for(size_t k = 0; k < K; k++){
		for(auto i: _K){
			if(k == i){
			 IsPresent = true;
				break;
			}
			else{
				IsPresent = false;
			}
		}

		if(!IsPresent){
			ClustersNew.insert({NewKey, Clusters[k]});
			NewKey++;
		}
	}

	Clusters.swap(ClustersNew);
    K = K - _K.size();
}

//Ritorna il numero di cluster attuale
template<unsigned int DIM> unsigned int CategoricalModel<DIM>::ViewK() const{
    return K;
}

//Id dei cluster presenti attualmente
template<unsigned int DIM> void CategoricalModel<DIM>::ViewKey(vector<unsigned int>& Key) const{

	for(auto it=Clusters.cbegin(); it!=Clusters.cend(); it++){
		Key.push_back(it->first);
	}

}

//Mostra il peso globale dei cluster presenti attualemente
template<unsigned int DIM> void CategoricalModel<DIM>::ViewBeta(vector<double>& _AllBeta) {

	 if( !(_AllBeta.empty())){
	    std::cerr<< " error in CategoricalModel.ViewBeta: _AllBeta must be empty (but memory reserved) "<<std::endl;
		exit(1);
	 }

	double Betak=0.0;

	for(size_t i=0; i < K; i++){
        Betak = Clusters[i].ViewBeta();
		_AllBeta.push_back(Betak);
	}
}

//Mostra il peso globale del sub-cluster sinistro dei clusters attualemente presenti
template<unsigned int DIM> void CategoricalModel<DIM>::ViewBetaLeft(vector<double>& _AllBetaLeft){

     if( !(_AllBetaLeft.empty())){
	    std::cerr<< " error in CategoricalModel.ViewBetaLeft: _AllBetaLeft must be empty (but memory reserved) "<<std::endl;
		exit(1);
	 }

	 double BetaLeftk=0.0;

	 for(size_t i=0; i < K; i++){
        BetaLeftk = Clusters[i].ViewBetaLeft();
		_AllBetaLeft.push_back(BetaLeftk);
	}
}

//Mostra il peso globale del sub-cluster destro dei clusters attualemente presenti
template<unsigned int DIM> void CategoricalModel<DIM>::ViewBetaRight(vector<double>& _AllBetaRight){

     if( !(_AllBetaRight.empty())){
	    std::cerr<< " error in CategoricalModel.ViewBetaRight: _AllBetaRight must be empty (but memory reserved) "<<std::endl;
		exit(1);
	 }

	 double BetaRightk=0.0;

	 for(size_t i=0; i < K; i++){
        BetaRightk = Clusters[i].ViewBetaRight();
		_AllBetaRight.push_back(BetaRightk);
	}
}

//Stampa su file i parametri latenti dei cluster
template<unsigned int DIM>
void CategoricalModel<DIM>::PrintTheta(const std::string& FileName) {

    std::ofstream file2(FileName,std::ios::binary | std::ios::app );

    if (file2.fail()){
		std::cerr <<"ERROR GENERATED (LastTheta)!" <<std::endl;
		std::cerr <<"Cannot open LastTheta" <<std::endl;
		exit(1);
	}

	//std::ofstream file_theta("WriteTheta.txt",std::ios::app);

	THETA Theta;
	for(size_t k=0; k<K; k++){
        Clusters[k].ViewTheta(Theta);

		/*for(auto i: Theta)
			file_theta<<i<<" ";
		file_theta<<std::endl;
		*/

        file2.write(reinterpret_cast<char*>(Theta.data()),Theta.size()*sizeof(double));
        file2.flush();

	}

	//file_theta<<"################"<<std::endl;

    file2.close();
	//file_theta.close();
}

//Stampa a video le informazioni degli iperparametri dei parametri latenti
template<unsigned int DIM>
void CategoricalModel<DIM>::PrintLambdaInfo() const {
    std::cout<<"Lambda = "<<Lambda[0]<<std::endl;
    std::cout<<std::endl;

}

/*
template<unsigned int DIM>
void CategoricalModel<DIM>::ViewModel(){
	vector<double> T;
	vector<unsigned int> C;

	for(auto iter = Clusters.begin(); iter != Clusters.end(); iter ++){
		std::cout<<"CLUSTER : "<<iter->first<<std::endl;
		std::cout<<"Beta: "<<iter->second.ViewBeta();
		std::cout<<"BetaLeft: "<<iter->second.ViewBetaLeft();
		std::cout<<"BetaRight: "<<iter->second.ViewBetaRight();
		std::cout<<"Theta: "<<std::endl;
		iter->second.ViewTheta(T);
		for(auto i: T)
			std::cout<<i<<" ";
		std::cout<<std::endl;
		std::cout<<"ThetaLeft"<<std::endl;
		iter->second.ViewThetaLeft(T);
		for(auto i: T)
			std::cout<<i<<" ";
		std::cout<<std::endl;
		std::cout<<"ThetaRight"<<std::endl;
		iter->second.ViewThetaRight(T);
		for(auto i: T)
			std::cout<<i<<" ";
		std::cout<<std::endl;
		std::cout<<"WordCount"<<std::endl;
		iter->second.ViewStatistics(C);
		for(auto i: C)
			std::cout<<i<<" ";
		std::cout<<std::endl;
		std::cout<<"WordCountLeft"<<std::endl;
		iter->second.ViewStatisticsLeft(C);
		for(auto i: C)
			std::cout<<i<<" ";
		std::cout<<std::endl;
		std::cout<<"WordCountRight"<<std::endl;
		iter->second.ViewStatisticsRight(C);
		for(auto i: C)
			std::cout<<i<<" ";
		std::cout<<std::endl;
		std::cout<<"GlobalTable: "<<iter->second.ViewGlobalTable()<<std::endl;
		std::cout<<"GlobalTableLeft: "<<iter->second.ViewGlobalTableLeft()<<std::endl;
		std::cout<<"GlobalTableRight: "<<iter->second.ViewGlobalTableRight()<<std::endl;

		}
	}
*/

#endif // _MODEL_HPP_
