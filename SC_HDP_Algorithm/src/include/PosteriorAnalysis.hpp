#ifndef __POSTERIOR_ANALYSIS__
#define __POSTERIOR_ANALYSIS__

#include <RInside.h>
#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <math.h>
#include<tuple>
#include <iomanip>
#include "Type.hpp"

using std::deque;
using std::greater;
using std::vector;
using std::priority_queue;
using std::pair;
using std::unordered_map;
using std::tuple;


/*! \file PosteriorAnalysis.hpp
 *  Contains classes for the posterior analysis of the algorithm's results.
 *  For each model there is a derived and specialized class.
 *  The generic class provides the common interface to all model specific classes.
 *  \date February 2016
 */
 
 
namespace CoeffSimilitudine{

    template<typename T> double VectorNorm (vector<T>& v)
    {
        size_t dim = v.size();
        T norm = 0.0;
        for(size_t i=0; i<dim; i++)
            norm += v[i]*v[i];

        return std::sqrt(static_cast<double>(norm));
    }	
	
    template<typename T> T ScalarProdoct (vector<T>& v1 , vector<T>& v2)
    {
        size_t dim1 = v1.size();
        size_t dim2 = v2.size();

        if(dim1 != dim2){
            std::cerr<<"Error: I can't do a scalar prodoct of vectos with different dimensions"<<std::endl;
            exit(1);
        }

        T prod = 0.0;
        for(size_t i=0; i<dim1; i++)
            prod += v1[i]*v2[i];

        return prod;
    }
	
    template<typename T> double Coeff (vector<T>& v1 , vector<T>& v2)
    {
        double normv1 = VectorNorm(v1);
        double normv2 = VectorNorm(v2);

        T scalarv1v2 = ScalarProdoct(v1,v2);

        return(static_cast<double>(scalarv1v2)/(normv1*normv2));
    }

};

/*! \brief Sorting operator for <unsigned int,double> pairs.
*    The order is based on the second element in the pair.
*    \authors{Debora Parisi and Stefania Perego}
*	 \date February 2016
*/


struct greater_for_pair{

       bool operator() (const std::pair<unsigned int,double>& x,const std::pair<unsigned int,double>& y)const {return x.second > y.second;}
};

/*! \brief Generic class for the posterior analysis.
 *  Virtual class where all methods are null.
 *  Each inherited class must define all methods in the base class and, if necessary, add other methods.
 *  Calls R scripts.
 *  Computes LPML index, identifies the topics and detect the best clustering of data according to the least square criteria.
 *  \authors{Debora Parisi and Stefania Perego}
 *  \date February 2016
 */
template<typename Type, unsigned int DIM>
class GenericPosteriorAnalysis{

    public:
		
		/*! \brief Sets the MCMC chain of the number of clusters K
        */
		virtual void SetAllK() = 0;
		
		/*! \brief Sets the MCMC chain of the values for concentration parameter \f$ \alpha \f$
        */
		virtual void SetAllAlpha()= 0;
		
		/*! \brief Sets the MCMC chain of the values for concentration parameter \f$ \gamma \f$
        */
		virtual void SetAllGamma() = 0;
		
		/*! \brief Calls the R script for the analysis of the K chain
		 *  \param R - R istance
		 *  \param Burnin - number of initial values in the chain to be discarded
		 *  \param Thinning - keep a value in the chain every thinning values
		*/
		virtual void KPosteriorAnalysis(RInside&, const unsigned int, const unsigned int) = 0;
		
		/*! \brief Calls the R script for the analyis of Alpha and Gamma chains
		 *  \param R - R istance
		 *  \param AlphaBurnin - number of initial values in the Alpha chain to be discarded
		 *  \param AlphaThinning - keep a value in the Alpha chain every AlphaThinning values
		 *  \param GammaBurnin - number of initial values in the Gamma chain to be discarded
		 *  \param GammaThinning - keep a value in the Gamma chain every GammaThinning values
		 *  \param  AlphaTry - yes if you want to repeat the Alpha chain's analysis, no otherwise
		 *  \param  GammaTry - yes if you want to repeat the Gamma chain's analysis, no otherwise
		 */
		virtual void AGPosteriorAnalysis(RInside&, const unsigned int,const unsigned int,const unsigned int,
                                         const unsigned int, const char,const char) = 0;
		
		/*! \brief Sets the R working directory
		 * \param _wd - R working directory
		*/					
		virtual void Setwd(const std::string&) = 0;
		
		/*! \brief Sets the dimension of the hyperparameter fo the laten parameter's distribution
		 * \param _W - dimension of the hyperparameter fo the laten parameter's distribution
		*/  
		virtual void SetW(const unsigned int) = 0;
		
		/*! \brief Sets the number of groups
		*   \param _D - number of groups
		*/
		virtual void SetD(const unsigned int) = 0;
		
		/*! \brief Sets the total number of data
		*   \param _N - total number of data
		*/
		virtual void SetN(const unsigned int) = 0;
		
		/*! \brief Sets the number of iterations for the MCMC chain
		*	\param _Iterations - number of iterations
		*/
		virtual void SetIterations(const unsigned int) = 0;
		
		/*!	\brief Sets the number of iterations to discard in order to compute the LPML
		*	\param _burnin - iterations to discard in order to compute the LPML
		*/
		virtual void Setburnin(const unsigned int) = 0;
		
		/*! \brief Sets burnin and thinning for the K chain
		*	\param _Burnin - burnin
		*	\param _Thinning - thinning
		*/
		virtual void SetBT(const unsigned int, const unsigned) = 0;
		
		/*! \brief Finds the best clustering according to the least squares criteria
		*/ 
		virtual void LeastSquareClustering() = 0;

		/*! \brief Computes and prints on the terminal the LPML index
		*/ 
		virtual void LPML() = 0;
		
		/*! \brief Skims the Labels*.bin files and joins them in the Labels.bin file
		*	\param R - object of class RInside
		*/		
		virtual void UnioneLabels(RInside& R) = 0;
		
		/*! \brief Writes on binary files \f$ \theta \f$ and  \f$ \beta \f$ parameters of the optimal iteration
		*	\param file_nr - files in which looking for parameters
		*	\param move - in those files, number of initial iterations to skip
		*	\return number of clusters in the best clustering
		*/    
		virtual void WriteBestParams() = 0;
		
		/*! \brief Acquires documents' names
		*/   
		virtual void SetDocs() = 0;
		
		/*! \brief Associates documents to estimated topics
		*/    
		virtual void AssociatingDocs() = 0;

};

/*! \brief Class for the posterior analysis when data are categorical and \f$ H \f$ is the Dirichlet distribution.
 *  Reads and stores the results in suitable structures.
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */  

template <unsigned int DIM=1>
class CategoricalPosteriorAnalysis final: public GenericPosteriorAnalysis <TypeCategorical<DIM>,DIM>{

	public:
		
		/*! \brief Latent parameter: vector of topic's weights for the distinct words
         */
		using THETA = TypeCategorical<1>::THETA;   
		
		/*! \brief Topic id
        */
		using ClusterId = unsigned int;
	
		/*! \brief Document id
        */
		using GroupId = unsigned int;
	
		/*! \brief Word id
        */
		using DataId = unsigned int;
	
		/*! \brief Control variable
        */
		using Check = unsigned int;
		

	private:
	
		/*! \brief Corpus' vocabulary
        */
		unordered_map<DataId,std::string> Vocabulary;
		
		/*! \brief Containes the values taken on by \f$ \alpha \f$ in each iteration, when \f$ \alpha \f$ is random
        */
		vector<double> Alpha;   
		
		/*! \brief Containes the values taken on by \f$ \gamma \f$ in each iteration, when \f$ \gamma \f$ is random
        */
		vector<double> Gamma;
		
		/*! \brief Contains the number of topics inferred in each iteration
        */
		vector<unsigned int> AllK;                   
	
		/*! \brief R working directory
        */
		std::string wd;

		/*! \brief Number of distinct words
        */
		unsigned int W;
		
		/*! \brief Total number of words
        */
		unsigned int N;
		
		/*! \brief Number of documents
        */
		unsigned int D;
		
		/*! \brief Number of iteration in the MCMC chain
		*/	
		unsigned int Iterations;
		
		/*!	\brief Number of initial iteration to discard in order to compute the LPML
		*/
		unsigned int burnin;
		
		/*! \brief Number of initial iteration to discard in the K chain
		*/
		unsigned int Burnin;
		
		/*!	\brief Thinning for the K chain
		*/
		unsigned int Thinning;
		
		/*! \brief Letf iterations after skimming the chain
		*/
		unsigned int Iter_after;
		
		/*! \brief Optimal iteration according to the least square clustering method
		*/
		unsigned int iter_opt;
		
		/*! \brief Files where looking for best parameters
		*/
		unsigned int file_nr;
		
		/*! \brief Number of iterations to be skipped in those files    
		*/
		unsigned int move;
		
		/*!	\brief Optimal number of topics
		*/
		unsigned int K_opt;
		
		/*! \brief Structure that stores the documents' names.
		*/
		vector<std::string> docs;

	public:

		/*! \brief Default constructor
        */
		CategoricalPosteriorAnalysis() = default;
		
		/*! \brief Destructor
        */
		~CategoricalPosteriorAnalysis()=default;

		/*! \brief Sets the chain of the number of topics
        */
		void SetAllK();
				
		/*! \brief Sets the chain of the values for \f$ \alpha \f$
        */
		void SetAllAlpha();
		
		/*! \brief Sets the chain of the values for \f$ \gamma \f$
        */
		void SetAllGamma();
		
		/*! \brief Sets the vocabulary
        */
		void SetVocabulary();
		
		/*! \brief Calls the R script for the analysis of the K chain
		 *  \param R - R istance
		 *  \param Burnin - number of initial values in the chain to be discarded
		 *  \param Thinning - keep a value in the chain every thinning values
		*/
		void KPosteriorAnalysis(RInside& R, const unsigned int Burnin, const unsigned int Thinning);
		
		/*! \brief Calls the R script for the analyis of Alpha and Gamma chains
		 *  \param R - R istance
		 *  \param AlphaBurnin - number of initial values in the Alpha chain to be discarded
		 *  \param AlphaThinning - keep a value in the Alpha chain every AlphaThinning values
		 *  \param GammaBurnin - number of initial values in the Gamma chain to be discarded
		 *  \param GammaThinning - keep a value in the Gamma chain every GammaThinning values
		 *  \param  AlphaTry - yes if you want to repeat the Alpha chain's analysis, no otherwise
		 *  \param  GammaTry - yes if you want to repeat the Gamma chain's analysis, no otherwise
		 */
		void AGPosteriorAnalysis(RInside& R, const unsigned int AlphaBurnin,const unsigned int AlphaThinning,const unsigned int GammaBurnin,
                                 const unsigned int GammaThinning, const char AlphaTry,const char GammaTry);
		
		/*! \brief Sets the R working directory
		 * \param _wd - R working directory
		*/
		void Setwd(const std::string& _wd);
		
		/*! \brief Sets the number of distintc words
		 * \param _W - number of distintc words
		*/
		void SetW(const unsigned int _W);
		
		/*! \brief Sets the number of documents
		*   \param _D - number of documents
		*/
		void SetD(const unsigned int _D);
		
		/*! \brief Set the total number of words
		*   \param _N - total number of words
		*/
		void SetN(const unsigned int _N);
		
		/*! \brief Skims the Labels*.bin files and joins them in the Labels.bin file
		*	\param R - object of class RInside
		*/			
		void UnioneLabels(RInside& R);
		
		/*! \brief Computes least-square clustering
		*/
		void LeastSquareClustering();
		
		/*! \brief Writes on binary files \f$ \theta \f$ and  \f$ \beta \f$ parameters of the optimal iteration
		*	\param file_nr - files in which looking for parameters
		*	\param move - in those files, number of initial iterations to skip
		*	\return number of topics in the best clustering
		*/ 
		void WriteBestParams();
		
		/*! \brief Draws the wordclouds representing the estimated topics in the optimal iteration
		* 	\param R - object of class RInside
		*/		
		void WordClouds(RInside& R);
		
		/*! \brief Acquires documents' names
		*/
		void SetDocs();
		
		/*! \brief Associates documents to estimated topics
		*/
		void AssociatingDocs();
		
		/*! \brief Sets the number of iterations for the MCMC chain
		*	\param _Iterations - number of iterations
		*/	
		void SetIterations(const unsigned int _Iterations);
		
		/*!	\brief Sets the number of iterations to discard in order to compute the LPML
		*	\param _burnin - iterations to discard in order to compute the LPML
		*/
		void Setburnin(const unsigned int _burnin);
		
		/*! \brief Sets burnin and thinning for the K chain
		*	\param _Burnin - burnin
		*	\param _Thinning - thinning
		*/
		void SetBT(const unsigned int _Burnin, const unsigned int _Thinning);

		/*! \brief Computes and prints on the terminal the LPML index
		*/ 
		void LPML();
		
};

//Imposta la working directory di R
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::Setwd(const std::string& _wd){
    wd = _wd;
}

//Imposta la dimensione dell'iperparametro della distribuzione del  parametro latente
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetW(const unsigned int _W){
    W = _W;
}

//Imposta il numero di documenti
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetD(const unsigned int _D){
    D = _D;
}

//Imposta il numero totale di dati
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetN(const unsigned int _N){
    N = _N;
}

//Imposta la catena AllK
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetAllK(){

    std::ifstream file("./cpp_results/AllK1.txt");

	if (file.fail()){
		std::cerr <<"Error in SetAllK" <<std::endl;
		std::cerr <<"I cannot open AllK.txt" <<std::endl;
		exit(1);
	}

	std::string ss;
	unsigned int temp;

	while(std::getline(file,ss))
	{
        std::istringstream SSTR(ss);
        SSTR>>temp;
        AllK.push_back(temp);
    }
	
	if(AllK.size() != Iterations){
		std::cerr<<" number of iteration in the K chain is not equal to Iterations "<<std::endl;
		exit(1);
	}


}

//Imposta la catena Alpha
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetAllAlpha(){

    std::ifstream file ("./cpp_results/AllAlpha1.txt");

	if (file.fail()){
		std::cerr <<"Error in SetAllAlpha" <<std::endl;
		std::cerr <<"I cannot open AllAlpha.txt" <<std::endl;
		exit(1);
	}

	std::string ss;
    double temp;

	while(std::getline(file,ss)){
        std::istringstream SSTR(ss);
        SSTR>>temp;
        Alpha.push_back(temp);
    }
	
	if(Alpha.size() != Iterations){
		std::cerr<<" number of iteration in the alpha chain is not equal to Iterations "<<std::endl;
		exit(1);
	}

}

//Imposta la catena AllK
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetAllGamma(){

    std::ifstream file ("./cpp_results/AllGamma1.txt");

	if (file.fail()){
		std::cerr <<"Error in SetAllGamma" <<std::endl;
		std::cerr <<"I cannot open AllGamma.txt" <<std::endl;
		exit(1);
	}

	std::string ss;
    double temp;

	while(std::getline(file,ss)){
        std::istringstream SSTR(ss);
        SSTR>>temp;
        Gamma.push_back(temp);
    }
	
	if(Gamma.size() != Iterations){
		std::cerr<<" number of iteration in the gamma chain is not equal to Iterations "<<std::endl;
		exit(1);
	}


}


//Imposta il vocabolario
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetVocabulary(){

    std::ifstream file("../Vocabulary.txt");

    if(file.fail()){
        std::cerr<<"Error in SetVocabulary"<<std::endl;
        std::cerr<<"I cannot open Vocabulary.txt"<<std::endl;
        exit(1);
    }

	std::string ss,ss1;
	DataId id = 1;
	std::string word;

	while(std::getline(file,ss)){
		std::istringstream SSTR(ss);
        while(std::getline(SSTR,ss1,' ')){
            std::istringstream SSTR1(ss1);
            SSTR1 >> word;
            Vocabulary.insert({id,word});
            ++id;
		}
	}

}

//Chiama lo script R per l'analisi della catena dei K
template<unsigned int DIM>
void CategoricalPosteriorAnalysis<DIM>::KPosteriorAnalysis(RInside& R, const unsigned int Burnin, const unsigned int Thinning){

	try
	{

		 //Devo andare nella cartella base per poter caricare il sorgente RAnalysis.cpp
		 // Devo usare wd, variabile globale
    R.parseEvalQ(wd);
	//I'm loading the coda library without report the message when we load a library in R
	std::string str_lib = "suppressMessages(library(coda))";
	R.parseEvalQ(str_lib);

	unsigned int n = AllK.size();
	unsigned int iter =0;

        for(size_t i = Burnin; i<n; i=i + Thinning)
		iter++;

	Rcpp::NumericVector R_K(iter);

	iter = 0;
	for(unsigned int i = Burnin; i<n; i=i + Thinning){
		R_K[iter] = AllK[i];
		iter++;
	}

	R["N"] = n;
	R["after_burnin"] = iter;
	R["burnin"] = Burnin;
	R["thinning"] = Thinning;
	R["AllK"] = R_K;

        // Chiamo lo script
        std::string src = "source(\"../../src/Rscript/KPosteriorAnalysis.R\")";   // ./ vuol dire "dove sono ora" cioè nella cartella ProvaR
        R.parseEval(src);
	}
	catch(std::exception& ex)
	{
		std::cerr << "Exception caught: " << ex.what() << std::endl;

	}
	catch(...)
	{
		std::cerr << "Unknown exception caught" << std::endl;
	}
	return;
}

//Chiama lo script R per l'analisi delle catene Alpha e Gamma
template<unsigned int DIM>
void CategoricalPosteriorAnalysis<DIM>::AGPosteriorAnalysis(RInside& R, const unsigned int AlphaBurnin,const unsigned int AlphaThinning,
                                                   const unsigned int GammaBurnin, const unsigned int GammaThinning, const char AlphaTry,const char GammaTry){

	try
	{

    R.parseEvalQ(wd);
	std::string str_lib = "suppressMessages(library(coda))";
	R.parseEvalQ(str_lib);

	R["alpha_try"] = AlphaTry;
	R["gamma_try"] = GammaTry;

	if(AlphaTry == 'y'){

		unsigned int n_Alpha = Alpha.size();
		unsigned int iter_Alpha =0;

        	for(size_t i = AlphaBurnin; i<n_Alpha; i=i + AlphaThinning)
			iter_Alpha++;

		Rcpp::NumericVector R_Alpha(iter_Alpha);

		iter_Alpha = 0;
		for(size_t i = AlphaBurnin; i<n_Alpha; i=i + AlphaThinning){
			R_Alpha[iter_Alpha] = Alpha[i];
			iter_Alpha++;
		}

		R["N_alpha"] = n_Alpha;
		R["after_alphaburnin"] = iter_Alpha;
		R["alphaburnin"] = AlphaBurnin;
		R["alphathinning"] = AlphaThinning;
		R["AllAlpha"] = R_Alpha;
	}
	if(GammaTry == 'y'){

		unsigned int n_Gamma = Gamma.size();
		unsigned int iter_Gamma =0;

        	for(unsigned int i = GammaBurnin; i<n_Gamma; i=i + GammaThinning)
			iter_Gamma++;

		Rcpp::NumericVector R_Gamma(iter_Gamma);

		iter_Gamma = 0;
		for(size_t i = GammaBurnin; i<n_Gamma; i=i + GammaThinning){
			R_Gamma[iter_Gamma] = Gamma[i];
			iter_Gamma++;
		}

		R["N_gamma"] = n_Gamma;
		R["after_gammaburnin"] = iter_Gamma;
		R["gammaburnin"] = GammaBurnin;
		R["gammathinning"] = GammaThinning;
		R["AllGamma"] = R_Gamma;
	}

        // Chiamo lo script
        std::string src = "source(\"../../src/Rscript/AGPosteriorAnalysis.R\")"; 
        R.parseEval(src);
	}
	catch(std::exception& ex)
	{
		std::cerr << "Exception caught: " << ex.what() << std::endl;

	}
	catch(...)
	{
		std::cerr << "Unknown exception caught" << std::endl;
	}

}



//Calcola e stampa a video l'indice LPML
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::LPML(){

    unsigned int MaxIt = Iterations - burnin ;
	
	std::ifstream lpml ("../cpp_results/CPO.bin");
    if (lpml.fail()){
		std::cerr <<"Error in LPML" <<std::endl;
		std::cerr <<"I cannot open CPO" <<std::endl;
		exit(1);
	}

	double LPML = 0.0;

    double *tmp = new double [N];
    lpml.read(reinterpret_cast<char*>(tmp), N*sizeof(double));
    vector<double> CPO (tmp, tmp + N);
    delete [] tmp;

	for_each(CPO.begin(), CPO.end(),[](double &den){den = 1.0/den;});
	for_each(CPO.begin(), CPO.end(), [MaxIt](double &den){ den = 1.0/((1.0/MaxIt)*den); });

    for (vector<double>::size_type i=0; i<CPO.size(); i++)
        LPML += std::log(CPO[i]);

    std::cout<<"LPML = "<<LPML<<std::endl;
}


// unisce le etichette
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::UnioneLabels(RInside& R){
	
    std::string FileNameLabels = "../cpp_results/Labels";
    std::string Exe = ".bin";
    std::string FileLables;

    std::streamsize howmove_labels;	
    unsigned int Iter = Iterations - 500;
    unsigned int counter = Thinning -1; // ogni thinning deve riazzerarsi
	unsigned int iter_after = 0; //numero iterazioni rimaste dopo la scrematura
		
    std::ofstream Labels ("../Labels.bin", std::ios::binary);
   
    if(Labels.fail()){

        std::cerr<<"I can't open Labels.bin"<<std::endl;
        exit(1);

    }
   
    for(size_t it=Burnin; it<=Iter; it=it + 500){

        std::stringstream add;

        add<<FileNameLabels<<it<<Exe;
        FileLables.clear();
        FileLables = add.str();

        std::cout<<"Nome del file "<<FileLables<<std::endl;

        std::ifstream templabel  (FileLables,std::ios::binary);


        if(templabel.fail()){
            std::cerr<<"Error. I can't open the minifiles"<<FileLables<<std::endl;
            exit(1);
        }


        for(size_t k=it; k<it+500; k++){

            counter++;


                if(Thinning != counter){
                    howmove_labels = sizeof(unsigned int)*N;
                    templabel.seekg(howmove_labels,std::ios::cur);
                }
                else{
                    unsigned int *tmp_labels = new unsigned int [N];
                    templabel.read(reinterpret_cast<char*>(tmp_labels), sizeof(unsigned int)*N);
                    std::vector<unsigned int> tmpL (tmp_labels,tmp_labels + N);


                Labels.write(reinterpret_cast<char*>(tmpL.data()), sizeof(unsigned int)*N);
				++ iter_after;
                }
                if(Thinning == counter){
                    counter = 0;
                }
    
        }

        templabel.close();
    }

	Labels.close();
		
	Iter_after = iter_after;
	
}




template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::LeastSquareClustering(){


	vector<vector<unsigned int>> Labels;
	
	//// lettura etichette ////
	
	std::ifstream file("../Labels.bin",std::ios::binary);

	if(file.fail()){
		std::cerr<<"Error in LoadLabels: cannot open Labels.bin"<<std::endl;
		exit(1);
	}
	
	vector<unsigned int> temp_label;
	
	//file.seekg(BurnIn*N*sizeof(unsigned int),std::ios::beg);
	
	for(size_t it = 0; it < Iter_after; ++it){

		unsigned int* tmp = new unsigned int[N];

		file.read(reinterpret_cast<char*>(tmp), N*sizeof(unsigned int));

		temp_label.assign(tmp,tmp + N);

		Labels.push_back(temp_label);

		delete [] tmp;
		
		//if(Thinning>1) file.seekg((Thinning-1)*N*sizeof(unsigned int),std::ios::cur);

	}
	
	
	file.close();
    
    //vettore in cui vengono memorizzate le somme;
    
    vector<double> ss(Labels.size(),0.0);
	
	//std::cout<<"metodo"<<std::endl;

        
    // calcolo le somme per ogni iterazione
    
    double pp;
    unsigned int count = 0;
    double deviance;
	size_t i,j,iter;
	unsigned int N_par = N;
    
    //double start = omp_get_wtime();

    #pragma omp parallel default(none) private(i,j,iter,pp,deviance) firstprivate(count) shared(ss,N_par,Labels,std::cout)
    {
		vector<double> ss_private(Labels.size(),0.0);
    	
		#pragma omp for schedule(dynamic,1)
		for(i = 0; i<N_par; ++i){

			for(j=0; j<i; ++j){
				
				//calcolo pp[i,j]
                for(iter=0; iter<Labels.size(); ++iter){
                    
                    if(Labels[iter][i]==Labels[iter][j]) ++count;
                }
				pp = count / static_cast<double>(Labels.size());
                
				
				for(iter=0; iter<Labels.size(); ++iter){
					deviance = (Labels[iter][i] == Labels[iter][j] ? 1 : 0 ) - pp;
					ss_private[iter] += deviance*deviance;			
				}
				
				count = 0;
							
			}
    	}
		
		#pragma omp critical
		{
			for(iter=0; iter<Labels.size(); ++iter)
				ss[iter] += ss_private[iter];
		}
	}

    
   // double stop = omp_get_wtime();
    
    //std::cout<<"time "<<stop-start<<std::endl;
	
	
	std::ofstream file_out("../LeastSquare.txt");
	
	
	auto iter_value = std::min_element(ss.cbegin(),ss.cend());
    double value = *iter_value;
    unsigned int pos = std::distance(ss.cbegin(),iter_value);
    
    file_out<<"value "<<value<<" pos "<<pos<<std::endl;
    
    for(size_t it=0; it<Labels.size();++it)
        if((it!=pos) && (ss[it]==value))
            file_out<<"another position found "<<it<<std::endl;
    
	
    unsigned int K = *(std::max_element(Labels[pos].cbegin(), Labels[pos].cend()))+1;

	vector<unsigned int> BestClustering(K,0);

	for(std::vector<unsigned int>::const_iterator it = Labels[pos].cbegin(); it != Labels[pos].cend(); ++it)
		++BestClustering[*it];
    
    file_out<<std::endl;
	file_out<<"Best clustering found at iteration "<<pos<<std::endl;
	file_out<<std::endl;
	file_out<<"Best clustering has "<<K<<" clusters"<<std::endl;
	file_out<<std::endl;
	for(unsigned int k = 0; k<K;++k)
		file_out<<"Cluster "<<k<<" has "<<BestClustering[k]<<" elements"<<std::endl;
	file_out<<std::endl;
	
	iter_opt = pos;
	
	unsigned int iter_corrisp = (iter_opt * Thinning) + Burnin;
	unsigned int migl = (iter_corrisp / 1000 ) *1000; // 2000
	unsigned int cent = iter_corrisp - migl; // 93
	
	if(cent < 500){
		file_nr = migl;
		move = cent;
	}else{
		file_nr = migl +500;
		move = cent - 500;
		
	}

}



template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::WriteBestParams(){
	
	
	std::string FileNameTheta = "../cpp_results/Theta";
	std::string FileNameBeta = "../cpp_results/Beta";
	std::string FileNamePi = "../cpp_results/Pi";
    std::string Exe = ".bin";
    std::string FileTheta;
	std::string FileBeta;
	std::string FilePi;
	
	std::stringstream add_theta;
	
	add_theta<<FileNameTheta<<file_nr<<Exe;
    FileTheta = add_theta.str();

    std::cout<<"Nome del file "<<FileTheta<<std::endl;
	
	std::stringstream add_beta;
	
	add_beta<<FileNameBeta<<file_nr<<Exe;
    FileBeta = add_beta.str();

    std::cout<<"Nome del file "<<FileBeta<<std::endl;
	
	std::stringstream add_pi;
	
	add_pi<<FileNamePi<<file_nr<<Exe;
    FilePi = add_pi.str();

    std::cout<<"Nome del file "<<FilePi<<std::endl;

	
	std::ifstream file_theta(FileTheta, std::ios::binary);
	std::ifstream file_beta(FileBeta, std::ios::binary);
	std::ifstream file_pi(FilePi, std::ios::binary);
	std::ofstream out_theta("../BestTheta.bin", std::ios::binary);
	std::ofstream out_beta("../BestBeta.bin", std::ios::binary);
	std::ofstream out_pi("../BestPi.bin", std::ios::binary);
	
	if(file_theta.fail()){
		std::cerr<<"cannot open Theta"<<std::endl;
		exit(1);
	}

	if(file_beta.fail()){
		std::cerr<<"cannot open Beta"<<std::endl;
		exit(1);
	}
	
	if(file_pi.fail()){
		std::cerr<<"cannot open Pi"<<std::endl;
		exit(1);
	}
	
	if(out_beta.fail()){
		std::cerr<<"cannot open BestBeta"<<std::endl;
		exit(1);
	}
	
	if(out_theta.fail()){
		std::cerr<<"cannot open BestTheta"<<std::endl;
		exit(1);
	}
	
	if(out_pi.fail()){
		std::cerr<<"cannot open BestPi"<<std::endl;
		exit(1);
	}
	
	
	 unsigned int K = AllK[file_nr + move];
	
	//std::cout<<"K "<<K<<std::endl;
	
	unsigned int cum_sum = 0;
	
	for(size_t i=file_nr; i<file_nr+move; ++i)
		cum_sum += AllK[i];
	
	//std::cout<<"cum "<<cum_sum<<std::endl;
	
	std::streamsize how_move_beta = cum_sum * sizeof(double);
	
	file_beta.seekg(how_move_beta,std::ios::beg);
	
	vector<double> beta;
	
	double* tmp_beta = new double[K];
	
	file_beta.read(reinterpret_cast<char*>(tmp_beta), K*sizeof(double));
	
	beta.assign(tmp_beta,tmp_beta + K);
	
	out_beta.write(reinterpret_cast<char*>(beta.data()), sizeof(double)*K);
	
	beta.clear();
	
	delete [] tmp_beta;
	
	//unsigned int W = 100;
	
	std::streamsize how_move_theta = cum_sum * W * sizeof(double);
	
	file_theta.seekg(how_move_theta,std::ios::beg);
	
	vector<double> theta;
	
	for(size_t k=0; k<K; ++k){
		double* tmp_theta = new double[W];
		file_theta.read(reinterpret_cast<char*>(tmp_theta), W*sizeof(double));
		theta.assign(tmp_theta, tmp_theta + W);
		out_theta.write(reinterpret_cast<char*>(theta.data()), sizeof(double)*W);
		delete [] tmp_theta;
		theta.clear();
			
	}
	
	std::streamsize how_move_pi = cum_sum * D * sizeof(double);
	
	file_pi.seekg(how_move_pi,std::ios::beg);
	
	vector<double> pi;
	
	for(size_t d=0; d<D; ++d){
		double* tmp_pi = new double[K];
		file_pi.read(reinterpret_cast<char*>(tmp_pi), K*sizeof(double));
		pi.assign(tmp_pi, tmp_pi + W);
		out_pi.write(reinterpret_cast<char*>(pi.data()), sizeof(double)*K);
		delete [] tmp_pi;
		pi.clear();
			
	}	

    file_theta.close();
	file_beta.close();
	file_pi.close();
	out_theta.close();
	out_beta.close();
	out_pi.close();
	
	K_opt = K;
	
}


template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::WordClouds(RInside& R){
	
	vector<pair<unsigned int,double>> Beta_opt;
		
	vector<pair<unsigned int,vector<double>>> Theta_opt;
			
	std::ifstream LastBeta("../BestBeta.bin",std::ios::binary);

    if(LastBeta.fail()){
		std::cerr<<"Doesn't open file: BestBeta"<<std::endl;
		exit(1);
    }
	
	
	double *tmp_beta = new double [K_opt];
    LastBeta.read(reinterpret_cast<char*>(tmp_beta), K_opt*sizeof(double));
    std::vector<double> TempBeta(tmp_beta, tmp_beta + K_opt);
    delete [] tmp_beta;
	
	for(size_t k=0; k<K_opt; ++k)
		Beta_opt.push_back({k,TempBeta[k]});
	
	
	std::ifstream LastTheta("../BestTheta.bin",std::ios::binary);

    if(LastTheta.fail()){
		std::cerr<<"Doesn't open file: BestTheta"<<std::endl;
		exit(1);
    }

	
	for(size_t k=0; k<K_opt; ++k){
        double *tmp_theta=new double[W];
        LastTheta.read(reinterpret_cast<char*>(tmp_theta),W*sizeof(double));
        std::vector<double> fun(tmp_theta, tmp_theta + W);
        delete [] tmp_theta;
        Theta_opt.push_back({k,fun});
    }

    char Try='y'; 
    unsigned int w; 
	unsigned int n;

    std::cout<<std::endl;
    std::cout<<"**********************************************"<<std::endl;
    std::cout<<"             Creating WordClouds               "<<std::endl;
    std::cout<<"**********************************************"<<std::endl;
    std::cout<<std::endl;
    while(Try == 'y'){
        std::cout<<std::endl;
        std::cout<<"How many topics do you want to visualize?  Insert a number between 1 - "<<K_opt<<": ";
        std::cin>>n;

        while(n>K_opt){
            std::cout<<"insert again a number between 1 - "<<K_opt<<": ";
            std::cin>>n;
        }
        std::cout<<std::endl;
        std::cout<<"How many words in each topics do you want to visualize (Min 3, Max "<<W<<")?  w = "<<std::endl;
        std::cin>>w;

        while(w>W){
            std::cout<<"insert again a number between 1 - "<<W<<": ";
            std::cin>>w;
        }
        std::cout<<std::endl;
        
		
		//Riordino il vettore Beta_opt per sapere quali sono gli n topic che devo prendere
		
		sort(Beta_opt.begin(), Beta_opt.end(), greater_for_pair());
	
		
		//etichette dei topic da mostrare
		
		vector<unsigned int> labels;
		
		for(size_t i=0; i<n; ++i)
			labels.push_back(Beta_opt[i].first);
		
		
		//ora vado a prendere in Theta_opt questi n topic, devo riordinare i pesi in senso decrescente e prendere le prime w parole
		//ma poichè mi serve tenere traccia degli id delle parole che prendo, mi serve un contenitore ausiliario, coppie id-peso parola
		
		vector<pair<unsigned int, double>> Theta_opt_aux;
		vector<vector<double>> Theta_opt_R;
		vector<double> temp_theta;
		vector<vector<unsigned int>> word_id;
		vector<unsigned int> temp_word_id;
		
		
		for(size_t i=0; i<n; ++i){
			
			temp_theta = Theta_opt[labels[i]].second;


			for(size_t id=0; id<W; ++id){
				Theta_opt_aux.push_back({id,temp_theta[id]});
			}
			
			sort(Theta_opt_aux.begin(),Theta_opt_aux.end(),greater_for_pair());

			temp_theta.clear();
			
			for(size_t id=0; id<w; ++id){
				temp_theta.push_back(Theta_opt_aux[id].second);
				temp_word_id.push_back(Theta_opt_aux[id].first+1);   // gli id nel vocabolario vanno da 1 a W
			}
			
			Theta_opt_R.push_back(temp_theta);
			word_id.push_back(temp_word_id);
		
			
			temp_word_id.clear();
			Theta_opt_aux.clear();
			
		}
		
        try {
		//Devo andare nella cartella base
        R.parseEvalQ(wd);
        std::string str_lib = "suppressMessages(library(wordcloud))";
        std::string str_lib2 = "suppressMessages(library(RColorBrewer))";
        R.parseEvalQ(str_lib);
        R.parseEvalQ(str_lib2);
        // Trasformo Beta dalla codifica in c++ a NumericVector, importabile in R
        //const unsigned int K = AllBeta.size();
		Rcpp::List RWeights(n);
		Rcpp::List RWords(n);
		Rcpp::CharacterVector RVocab(W);
		Rcpp::NumericVector RTopics(n);

		for(size_t k= 0; k<n; ++k)
			RWeights[k]= Theta_opt_R[k];

		
        for(size_t k= 0; k<n; ++k)
			RWords[k]= word_id[k];

        for(size_t It_idx = 0; It_idx<n; It_idx++)
            RTopics[It_idx] = labels[It_idx] + 1;

        for(unsigned int i=1; i<=W; i++)
            RVocab[i-1] = Vocabulary[i];

		// Ora posso mandare il vettore nel workspace di R

        R["Weights"] = RWeights;
        R["Dim"] = n;
        R["Words"] = RWords;
        R["Labels"] = RTopics;
        R["Dim_Voc"] = W;
        R["Voc"] = RVocab;
        R["Dim_w"] = w;


        // Chiamo lo script
			
		std::string src = "source(\"../../src/Rscript/WordCloud.R\")";	
        R.parseEval(src);
        }
        catch(std::exception& ex)
        {
            std::cerr << "Exception caught: " << ex.what() << std::endl;
        }
        catch(...)
        {
            std::cerr << "Unknown exception caught" << std::endl;
        }

        std::cout<<"Would you like to change something? press y [yes] or n [no] ";
        std::cin>>Try;

        while(Try !='y' && Try != 'n'){
            std::cout<<"press y [yes] or n [no] "<<std::endl;
            std::cin>>Try;
        }

        }
		
}


template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetDocs(){
	
		std::ifstream file_docs("../DocNames.txt");
	
		if(file_docs.fail()){
			std::cerr<<"Cannot open file: DocNames.txt"<<std::endl;
			exit(1);
   		 }
	
		std::string ss_docs;
	
		while(std::getline(file_docs,ss_docs))
				docs.push_back(ss_docs);
	
		if(docs.size() != D){
			std::cerr<<" read nr of documents different from D"<<std::endl;
			exit(1);
    	}

}


template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::AssociatingDocs(){
	
// acquisisco i Best Pi
	
	vector<vector<double>> Pi_opt;
	vector<double> Pi_opt_aux;
	std::ifstream file_pi("../BestPi.bin",std::ios::binary);
	
	for(size_t d=0; d<D; ++d){
		
		double* tmp_pi = new double[K_opt];
		file_pi.read(reinterpret_cast<char*>(tmp_pi), K_opt*sizeof(double));
		Pi_opt_aux.assign(tmp_pi, tmp_pi + K_opt);
		Pi_opt.push_back(Pi_opt_aux);
		delete [] tmp_pi;
		
	}
	
	
	
	//############ versione max pi ###########
	
	// per ogni j, cerco in \pi_j il peso massimo \pi_jk e associo al doc j il topic k-1
	// nel vettore corrisp, lungo K, l'elemento k è un vettore di doc associati al topic k,
	// se il vettore è vuoto vuol dire che a quel topic non ho associato nessun doc (controllo!)
	// in questo modo ogni doc è stato associato ad uno e un solo topic
	// alla fine per ogni topic k stampo i nomi dei doc che gli sono stati associati
	
	/*vector<vector<std::string>> topic_to_docs(K_opt);
	
	for(size_t d = 0; d< D; ++d){
		auto max_pi = std::max_element(Pi_opt[d].cbegin(),Pi_opt[d].cend());
		unsigned int position = std::distance(Pi_opt[d].cbegin(),max_pi);
		if(position != K_opt) //vuole dire che il peso più grande è quello del cluster vuoto, non dovrebbe succedere
			topic_to_docs[position].push_back(docs[d]);
		else{
			std::cerr<<"position = K_opt"<<std::endl;
			exit(1);
		}
	}
	
	
	std::ofstream file("DocAnalysis.txt");
	
	for(size_t k= 0; k< K_opt ; ++k){
		
		file << "########### Documents associated to topic "<<k+1<< " ###########"<<std::endl;
		file <<std::endl;
		
		if (topic_to_docs[k].empty()) {// se al topic k non sono stati associati documenti
			file<<"No document associated to topic "<<k+1<<std::endl;
			file<<std::endl;
	    }else{
			for(auto i : topic_to_docs[k] )
				file<<i<<std::endl;
			file<<std::endl;
			file<<std::endl;
		}
		
	}
	// ######################################################
		*/
	//################## versione soglia #########################
	
	// per ogni topic k, controllo per ogni doc j se \pi_jk è maggiore di una certa soglia
	// per il momento fisso la soglia, poi verrà decisa dall'utente e potrebbe variare per ogni topic (interazione utente)
	// se la condizione della soglia è vera, inserisco l'id del documento nel vettore topics_to_docs, poi è tutto come prima
	// con questo metodo non ho la certezza di associare ogni doc a qualche topic, perciò metto un flag che diventa vero se associo un doc ad almeno un topic


	//double threshold = 0.05;//1.0/K_opt;
	
	double threshold;
	
	std::cout<<"For each topic k in each document j I check if pi_jk is greater than a threshold."<<std::endl;
	std::cout<<"If so, I associate document j to topic k"<<std::endl;
	std::cout<<"Set a threshold for pi_jk:"<<std::endl;
	std::cin>>threshold;
	std::cout<<std::endl;
	
	vector<bool> flags(D);

	vector<vector<std::string>> topics_to_docs(K_opt);
	
	for(size_t k=0; k<K_opt; ++k){
		
		for(size_t d=0; d< D; ++d){
			
			if(Pi_opt[d][k] >= threshold){
				
				topics_to_docs[k].push_back(docs[d]);
				flags[d] = true;			
			}			
		}		
	}
	
	std::ofstream file("../DocAnalysis.txt");
	
	for(size_t k= 0; k< K_opt ; ++k){
		
		file << "########### Documents associated to topic "<<k+1<< " ###########"<<std::endl;
		file <<std::endl;
		
		if (topics_to_docs[k].empty()) {// se al topic k non sono stati associati documenti
			file<<"No document associated to topic "<<k+1<<std::endl;
			file<<std::endl;
	    }else{
			for(auto i : topics_to_docs[k] )
				file<<i<<std::endl;
			file<<std::endl;
		}
		
	}
	
	for(size_t d=0; d<D; ++d)
		if(flags[d] == false)
			std::cout<<"document "<<docs[d]<<" with id "<<d+1<<" not associated to any topic."<<std::endl;
	
	
	//###########################################################
	

}


// Imposta il numero di iterazioni
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetIterations(const unsigned int _Iterations){
	
	Iterations = _Iterations;

}


// Imposta il burnin per il calcolo di LPML
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::Setburnin(const unsigned int _burnin){
	
	burnin = _burnin;
}

// Imposta burnin e thinning per la catena dei K
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetBT(const unsigned int _Burnin, const unsigned int _Thinning){
	
	Burnin = _Burnin;
	Thinning = _Thinning;

}


#endif

