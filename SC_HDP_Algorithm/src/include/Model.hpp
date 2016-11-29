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
 * Class that manages all clusters, sampling equations and functions that depend on the chosen model.
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */

/*! \brief Interface for the Model class
 *  Abstract class where all methods are virtual.
 *  Classes that inherit from ModelGeneric sample latent parameters and manages related hyperparameters, which are specific to the chosen model.
 *  Manages model specific functions, such as likelihood, marginals and other densities.
 *  Removes and adds clusters.
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */
template<typename Type, unsigned int DIM=1>
class ModelGeneric{

    public:    

        /*! \brief Retrieves current number of clusters
         *  \return current number of clusters
         */
        virtual unsigned int ViewK() const = 0;

        /*! \brief Retrieves current clusters' ids
         *  \param Key - object where current clusters' ids are stored
         */
        virtual void ViewKey(vector<unsigned int>& ) const = 0;

        /*! \brief Fixes hyperparameters for the latent parameters' distribution
         *  \param _Lambda - hyperparameters for the latent parameters' distribution
         */
        virtual void SetHyperparameter(const typename Type::HYP&) = 0;

        /*! \brief Fixes hyperparameters for the latent parameters' distribution with default values
         *  \param W - hyperparameters' dimension
         */
        virtual void DefaultHyperparameter(size_t) = 0;

        /*! \brief Sets clusters assigning initial values for their weights
         *  \param _K - Initial number of clusters
         */
        virtual void SetInitialClusters(unsigned int) = 0;

        /*! \brief Computes the logarithm of the cluster's marginal likelihood, given the cluster id
         *  \param _K - cluster id
         */   
        virtual double Marginalized_Loglikelihood(const unsigned int)  = 0;

        /*! \brief Computes the loglikelihood of datum X, given the cluster id
         *  \param X - datum
         *  \param _K - cluster id
         *  \return loglikelihood of datum X in cluster k
         */
        virtual double Loglikelihood(const typename Type::Point , const unsigned int) = 0;

        /*! \brief Computes the loglikelihood of datum X, given the left subcluster id
         *  \param X - datum
         *  \param _K - left subcluster id
         *  \return loglikelihood of datum X in left subcluster of cluster k
         */
        virtual double LoglikelihoodLeft(const typename Type::Point, const unsigned int) = 0;

        /*! \brief Computes the loglikelihood of datum X, given the right subcluster id
         *  \param X - datum
         *  \param _K - right subcluster id
         *  \return loglikelihood of datum X in right subcluster of cluster k
         */
        virtual double LoglikelihoodRight(const typename Type::Point, const unsigned int) = 0;

        /*! \brief Computes the latent parameter's density, given the cluster id
         *  \param _K - cluster id
         *  \return logdensity of cluster k
         */
        virtual long double LogDensity(const unsigned int) = 0;

        /*! \brief Updates latent parameters of all clusters
         *  \param Gen - parallel random number generator
         */
        virtual void UpdateThetaCluster(omprng&) = 0;

        /*! \brief Updates latent parameters of left and right subclusters of all clusters
         *  \param Gen - parallel random number generator
         */
        virtual void UpdateThetaSubCluster(omprng&) = 0;

        /*! \brief Updates latent parameters of one cluster, given its id
         *  \param _K - cluster id
         *  \param Gen - parallel random number generator
         */
        virtual void UpdateOneThetaCluster(const unsigned int,omprng&) = 0;

        /*! \brief Updates latent parameters of left and right subclusters of one cluster, given its id
         *  \param _K - cluster id
         *  \param Gen - parallel random number generator
         */
        virtual void UpdateOneThetaSubCluster(const unsigned int,omprng&) = 0;

        /*! \brief Adds an empty cluster to the current clusters and updates the current number of clusters
         *  \param _K - new cluster id
         */  
        virtual void AddOneCluster(const unsigned int) = 0;

        /*! \brief Removes a cluster given its id and updates the current number of clusters
         *  \param _K - id of the cluster to be removed
         */  
        virtual void RemoveOneCluster (const unsigned int) = 0;

        /*! \brief Removes multiple clusters given their id and updates the current number of clusters
         *  \param _K - vector of clusters' ids to be removed
         */   
        virtual void RemoveClusters (const vector<unsigned int>&) = 0;

        /*! \brief Retrieves current clusters' global weights
         *  \param AllBeta - vector that will be filled with global weights
         */       
        virtual void ViewBeta(vector<double>&) = 0;

        /*! \brief Retrieves global weights for left subclusters of current clusters
         *  \param AllBetaLeft - vector that will be filled with global weights of left subclusters
         */   
        virtual void ViewBetaLeft(vector<double>&) = 0;

        /*! \brief Retrieves global weights for right subclusters of current clusters
         *  \param AllBetaRight - vector that will be filled with global weights of right subclusters
         */     
        virtual void ViewBetaRight(vector<double>&) = 0;

        /*! \brief Print to file values of current clusters' latent parameters
         */  
        virtual void PrintTheta(const std::string&) = 0;

        /*! \brief Print to screen information about hyperparameters of latent parameters' distribution
         */ 
        virtual void PrintLambdaInfo() const  = 0;

       // virtual void ViewModel() = 0;

};

/*! \brief Specialized class for topic modeling, where data are categorical and \f$ H \f$ is the Dirichlet distribution.
 *  This class is used to sample latent parameters from the Dirichlet distribution, to compute likelihood and marginal distribution for categorical data.
 *  It also removes and adds topics and print values of latent parameters on file.
 *  \authors{Debora Parisi and Stefania Perego}
 *  \date February 2016
 */
template <unsigned int DIM=1>
class CategoricalModel final: public ModelGeneric<TypeCategorical<DIM>,DIM>{

    public:

        //DEFIINIZIONE DEI TIPI PUBBLICI
        /*! \brief Latent parameter: weights of distintc words in a topic
         */
        using THETA = TypeCategorical<1>::THETA;

        /*! \brief A datum
         */
        using POINT = TypeCategorical<1>::Point;

        /*! \brief Statistics to updated hyperparameters for latent parameter's distribution, namely the number of data assigned to topic and subtopic
         */    
        using STAT  = TypeCategorical<1>::STAT;

        /*! \brief Hyperparameters for latent parameter's distribution
         */
        using HYP   = TypeCategorical<1>::HYP;

    private:

        //DEFINIZIONE DEI TIPI PRIVATI

        /*! \brief Vector of hyperparameters
         */
        HYP Lambda;

        /*! \brief Current number of clusters
         */
        unsigned int K;

        /*! \brief Container for clusters
         */
        unordered_map<unsigned int, CategoricalCluster> Clusters;

        /*! \brief Number of threads
         */
        unsigned int OMP_NUM_THREADS;

    public:

        //COSTRUTTORI E DISTRUTTORE

        /*! \brief Default constructor
         */
        CategoricalModel();

		/*! \brief Copy costructor
		 */
		CategoricalModel(const CategoricalModel& mod)= default;

		/*! \brief Move constructor
		 */
		CategoricalModel(const CategoricalModel&& mod);

		/*! \brief Destructor
		 */
        ~CategoricalModel()=default;

        //OPERATORI

        /*! \brief Operator to access directly to methods of cluster k
         *  \param _K - cluster id
         */
        CategoricalCluster& operator[] (unsigned int _K);

		/*! \brief Assignment operator
		 *  \param mod - object of class CategoricalModel
		 */
        CategoricalModel& operator=(const CategoricalModel &mod);

        /*! \brief Move assignement operator
		 *  \param mod - object of class CategoricalModel
		 */
        CategoricalModel& operator=(CategoricalModel &&mod);

        /*! \brief Fixes hyperparameters for the latent parameters' distribution
         *  \param _Lambda - hyperparameters for the latent parameters' distribution
         */
        void SetHyperparameter(const HYP& _Lambda);

        /*! \brief Fixes hyperparameters for the latent parameters' distribution with default values
         *  \param W - hyperparameters' dimension
         */
        void DefaultHyperparameter(size_t W);

        /*! \brief Sets topics assigning initial values for their weights
         *  \param _K - Initial number of topics
         */
		void SetInitialClusters(unsigned int _K);

        //METODI SPECIFICI DEL MODELLO:
        //VERISIMIGLIANZE MARGINALE DENSITÀ

       /*! \brief Computes the logarithm of the topic's marginal likelihood, given the topic id
         *  \param _K - topic id
         */   
		double Marginalized_Loglikelihood(const unsigned int _K);

        /*! \brief Computes the loglikelihood of datum X, given the topic id
         *  \param X - datum
         *  \param _K - topic id
         *  \return loglikelihood of datum X in topic k
         */
        double Loglikelihood(const POINT X,const unsigned int _K);

        /*! \brief Computes the loglikelihood of datum X, given the left subtopic id
         *  \param X - datum
         *  \param _K - left subtopic id
         *  \return loglikelihood of datum X in left subtopic of topic k
         */
        double LoglikelihoodLeft(const POINT X,const unsigned int _K);

        /*! \brief Computes the loglikelihood of datum X, given the right subtopic id
         *  \param X - datum
         *  \param _K - right subtopic id
         *  \return loglikelihood of datum X in right subtopic of topic k
         */
        double LoglikelihoodRight(const POINT X,const unsigned int _K);

        /*! \brief Computes the latent parameter's density, given the topic id
         *  \param _K - topic id
         *  \return logdensity of topic k
         */
        long double LogDensity(const unsigned int _K);

        //METODI PER AGGIORNARE I PARAMETRI LANTENTI

        /*! \brief Updates latent parameters of all topics
         *  \param Gen - parallel random number generator
         */
        void UpdateThetaCluster(omprng& Gen);

        /*! \brief Updates latent parameters of left and right subtopics of all topics
         *  \param Gen - parallel random number generator
         */
        void UpdateThetaSubCluster(omprng& Gen);

        /*! \brief Updates latent parameters of one topic, given its id
         *  \param _K - topic id
         *  \param Gen - parallel random number generator
         */
        void UpdateOneThetaCluster(const unsigned int _K,omprng& Gen);

        /*! \brief Updates latent parameters of left and right subtopics of one topic, given its id
         *  \param _K - topic id
         *  \param Gen - parallel random number generator
         */
        void UpdateOneThetaSubCluster(const unsigned int _K,omprng& Gen);

        //METODI PER AGGIUNGERE E RIMUOVERE I CLUSTER

        /*! \brief Adds an empty topic to the current topics and updates the current number of topics
         *  \param _K - new topic id
         */ 
        void AddOneCluster (const unsigned int _k);

        /*! \brief Removes a topic given its id and updates the current number of topics
         *  \param _K - id of the topic to be removed
         */
        void RemoveOneCluster(const unsigned int _K);

        /*! \brief Removes multiple topics given their id and updates the current number of topics
         *  \param _K - vector of topics' ids to be removed
         */ 
        void RemoveClusters(const vector<unsigned int>& _K);

        //ALTRI METODI
        /*! \brief Retrieves current number of topics
         *  \return current number of topics
         */
        unsigned int ViewK() const;

        /*! \brief Retrieves current topics' ids
         *  \param Key - object where current topics' ids are stored
         */
        void ViewKey(vector<unsigned int>& Key) const;

        /*! \brief Retrieves current topics' global weights
         *  \param AllBeta - vector that will be filled with global weights
         */ 
	   	void ViewBeta(vector<double>& _AllBeta);

        /*! \brief Retrieves global weights for left subtopics of current topics
         *  \param AllBetaLeft - vector that will be filled with global weights of left subtopics
         */
		void ViewBetaLeft(vector<double>& _AllBetaLeft) ;

        /*! \brief Retrieves global weights for right subtopics of current topics
         *  \param AllBetaRight - vector that will be filled with global weights of right subtopics
         */ 
		void ViewBetaRight(vector<double>& _AllBetaRight);

        /*! \brief Print to file values of current topics' latent parameters
         */
		void PrintTheta(const std::string& FileName) ;

        /*! \brief Print to screen information about hyperparameters of latent parameters' distribution
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
        std::cerr<<"In CategoricalModel can only be equal to 1"<<std::endl;
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
                        std::cerr<<"Lambda and vector of counts have different dimension"<<std::endl;
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
           	std::cerr<<"Error in CategoricalModel: UpdateThetaCluster"<<std::endl;
            std::cerr<<"Lambda and vector of counts have different dimension"<<std::endl;
            exit(1);
        }
        if(Lambda.size()!=cRight.size()){
           	std::cerr<<"Error in CategoricalModel: UpdateThetaCluster"<<std::endl;
            std::cerr<<"Lambda and vector of counts have different dimension"<<std::endl;
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
        std::cerr<<"Error in CategoricalModel: UpdateOneThetaCluster"<<std::endl;
        std::cerr<<"Lambda and vector of counts have different dimension"<<std::endl;
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
        std::cerr<<"Error in CategoricalModel: UpdateThetaCluster"<<std::endl;
        std::cerr<<"Lambda and vector of counts have different dimension"<<std::endl;
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
        std::cerr<<"Error in CategoricalModel: AddOneCluster"<<std::endl;
        std::cerr<<"Cannot generate a non-existing key"<<std::endl;
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
