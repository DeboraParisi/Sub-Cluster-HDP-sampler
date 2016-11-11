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
 *  Contiene le classi che gestiscono l'analisi a posteriori dei risultati prodotti dall'algoritmo.
 *  Per ogni modello esiste una classe derivata e specializzata.
 *  La classe generica fornisce l'interfaccia comune a tutte le classi specifiche di un modello.
 *  \date Febbraio 2016
 */
 
/*!  \brief Namespace per il calcolo del coefficiente di similitudine 
 *   Contiene tre funzioni template utili al calcolo del coefficiente di similitudine tra vettori di un qualsiasi tipo
 */

namespace CoeffSimilitudine{

	/*! \brief Calcola la norma di un vettore
	*  \param v - vettore di cui calcolare la norma
	*  \return Norma del vettore
	*/
    template<typename T> double VectorNorm (vector<T>& v)
    {
        size_t dim = v.size();
        T norm = 0.0;
        for(size_t i=0; i<dim; i++)
            norm += v[i]*v[i];

        return std::sqrt(static_cast<double>(norm));
    }
	
	
	/*!  \brief Calcola il prodotto scalare tra due vettori 
	 *   \param v1 - primo vettore
	 *   \param v2 - secondo vettore
	 *   \return Prodotto scalare tra i due vettori
	 */
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
	
	/*! \brief Calcola il coefficiente di similitudine tra due vettori 
	 *  \param v1 - primo vettore
	 *  \param v2 - secondo vettore
	 *  \return Coefficiente di similitudine
	 */
    template<typename T> double Coeff (vector<T>& v1 , vector<T>& v2)
    {
        double normv1 = VectorNorm(v1);
        double normv2 = VectorNorm(v2);

        T scalarv1v2 = ScalarProdoct(v1,v2);

        return(static_cast<double>(scalarv1v2)/(normv1*normv2));
    }

};

/*! \brief Classe generica per l'analisi a posteriori
 *  Classe virtuale dove tutti i metodi sono null.
 *  Ogni classe che eredita deve definire tutti i metodi della classe base e, se necessario, può aggiungere altri metodi.
 *  Invoca gli script R per l'analisi delle catene MCMC.
 *  Calcola l'indice LPML, riconosce i topic e individua il miglior clustering secondo il criterio dei minimi quadrati
 */

template<typename Type, unsigned int DIM>
class GenericPosteriorAnalysis{

    public:
		
		/*! \brief Imposta la catena AllK
		 *  \return Lunghezza della catena
        */
		virtual unsigned int SetAllK() = 0;
		
		/*! \brief Imposta la catena Alpha
		 *  \return Lunghezza della catena
        */
		virtual unsigned int SetAllAlpha()= 0;
		
		/*! \brief Imposta la catena Gamma
		 *  \return Lunghezza della catena
        */
		virtual unsigned int SetAllGamma() = 0;
		
		/*! \brief Chiama lo script R che visualizza le traiettorie della distribuzione dei cluster nel corpus
		 *  \param R - Istanza di R
		 *  \param BestClustering - iterazione a cui è stato individuato il least square clustering
        */
		//virtual void VisualizeBeta(RInside&,unsigned long)= 0;
		
		/*! \brief Chiama lo script R per l'analisi della catena dei K
		 *  \param R - Istanza di R
		 *  \param Burnin - numero di valori iniziali della catena catena da scartare
		 *  \param Thinning - tiene un valore della catena ogni thinning valori
		*/
		virtual void KPosteriorAnalysis(RInside&, const unsigned int, const unsigned int) = 0;
		
		/*! \brief Chiama lo script R per l'analisi delle catene Alpha e Gamma
		 *  \param R - Istanza di R
		 *  \param AlphaBurnin - numero di valori iniziali della catena Alpha da scartare
		 *  \param AlphaThinning - tiene un valore della catena Alpha ogni AlphaThinning valori
		 *  \param GammaBurnin - numero di valori iniziali della catena Gamma da scartare
		 *  \param GammaThinning - tiene un valore della catena Gamma ogni GammaThinning valori
		 *  \param  AlphaTry - yes se si vuole ripetere l'analisi della catena Alpha, no altrimenti
		 *  \param  GammaTry - yes se si vuole ripetere l'analisi della catena Gamma, no altrimenti
		 */
		virtual void AGPosteriorAnalysis(RInside&, const unsigned int,const unsigned int,const unsigned int,
                                         const unsigned int, const char,const char) = 0;
		
		/*! Imposta la working directory di R
		 * \param _wd - working directory di R
		*/					
		virtual void Setwd(const std::string&) = 0;
		
		/*! Imposta la dimensione dell'iperparametro della distribuzione del  parametro latente
		 * \param _W - dimensione dell'iperparametro della distribuzione del  parametro latente
		*/
		virtual void SetW(const unsigned int) = 0;
		
		/*! \brief Imposta il numero di gruppi
		*   \param _D - numero di gruppi
		*/
		virtual void SetD(const unsigned int) = 0;
		
		/*! \brief Imposta il numero totale di dati
		*   \param _N - numero totale di dati
		*/
		virtual void SetN(const unsigned int) = 0;
		
		/*! \brief Imposta la ricerca dei cluster, acquisendo da file i \f$ \theta \f$ e i  \f$ \beta \f$ delle ultime iterazioni
		 *  \param MaxIt - numero di iterazioni dell'algoritmo
		*/
		virtual void SetTheta(const unsigned long) = 0;
		
		/*! \brief Legge da file e imposta i \f$ \pi \f$
		 *  \param MaxIt - numero di iterazioni dell'algoritmo
		*/
		virtual void SetPi(const unsigned long) = 0;

		/*! \brief Acquisisce le etichette assegnate ai dati nelle ultime iterazioni
		 *  \param NrClusterings - numero di iterazioni da monitorare
		*/
		virtual void LoadLabels(const unsigned long) = 0;
		
		/*! \brief Individua il miglior clustering secondo il criterio dei minimi quadrati; si veda equazione (4.3) in Relazione_Parisi_Perego.pdf
		*   \return Iterazione a cui è stato individuato il miglior clustering
		*/
		virtual unsigned long LeastSquareClustering() = 0;

		/*! \brief Calcola e stampa a video l'indice LPML; si vedano equazioni (4.1)- (4.2) in Relazione_Parisi_Perego.pdf
		*   \param MaxIt - numero di iterazioni da considerare nel calcolo di LPML
		*/
		//virtual void LPML(unsigned long) = 0;
		
		/*! \brief Effettua il riconoscimento dei cluster
		*   \param R - Istanza di R
		*/
		//virtual void TrackingClusters(RInside&) = 0;

};

/*! Classe per l'analisi a posteriori specifica del modello Dirichlet-Categorical
 *  Legge e memorizza i risultati delle simulazioni in opportune strutture.
 */

template <unsigned int DIM=1>
class CategoricalPosteriorAnalysis final: public GenericPosteriorAnalysis <TypeCategorical<DIM>,DIM>{

	public:
		
		/*! \brief Parametro latente, vettore di pesi delle parole distinte nel cluster
         */
		using THETA = TypeCategorical<1>::THETA;
		
		/*! \brief Identificativo del cluster
        */
		using ClusterId = unsigned int;
	
		/*! \brief Identificativo del documento
        */
		using GroupId = unsigned int;
	
		/*! \brief Identificativo del dato
        */
		using DataId = unsigned int;
	
		/*! \brief Variabile di controllo
        */
		using Check = unsigned int;
		

	private:

		/*! \brief Struttura in cui memorizzare le informazioni sui topic riconosciuti.
		 *  La chiave della mappa è l'iterazione, il valore mappato è il blocco di topic individuati a quella iterazione.
		 *  Per ogni topic nel blocco si ha: l'etichetta assegnata al topic, il numero di parole in comune con i topic di altri blocchi aventi la stessa etichetta,
		 *  il peso globale del topic, il vettore di pesi delle parole che rappresenta il topic.
        */
		unordered_map<unsigned int,vector<tuple<GroupId,Check,double,vector<double>>>> Theta;
		
		/*! \brief Memorizza i vettori \f$ \beta \f$ delle ultime iterazioni
        */
		vector<vector<double>> AllBeta;
		
		/*! \brief Struttura in cui memorizzare i \f$ \pi \f$
		 * La chiave è l'iterazione, il valore mappato è un vettore lungo D, in posizione j ho il vettore dei pesi dei topic nel documento j
		*/
		unordered_map<unsigned int,vector<vector<double>>> AllPi;
		
		/*! \brief Struttura in cui memorizzo i Theta ottimali, con cui creare i wordcloud
		 * Nella posizione k del vettore ho il vettore \f$ \theta_k \f$
		*/
		vector<vector<double>> Theta_opt;
		
		/*! \brief Struttura in cui memorizzo i Pi ottimali, con cui associare i documenti ai topic
		 * Nella posizione j del vettore ho il vettore \f$ \pi_j \f$
		*/
		vector<vector<double>> Pi_opt;
		
		/*! \brief Struttura in cui memorizzo i Beta ottimali, con cui riordinare i topic dal più al meno importante
		 * Nella posizione k del vettore ho il peso \f$ \beta_k \f$
		*/
		vector<double> Beta_opt;
		
		/*! \brief Vocabolario del corpus
        */
		unordered_map<DataId,std::string> Vocabulary;
		
		/*! \brief Contiene i valori assunti da \f$ \alpha \f$ ad ogni iterazione, nel caso di prior su \f$ \alpha \f$
        */
		vector<double> Alpha;
		
		/*! \brief Contiene i valori assunti da \f$ \gamma \f$ ad ogni iterazione, nel caso di prior su \f$ \gamma \f$
        */
		vector<double> Gamma;
		
		/*! \brief Contiene il numero di cluster dedotti dall'algoritmo ad ogni iterazione
        */
		vector<unsigned int> AllK;
	
		/*! \brief Working directory per R
        */
		std::string wd;

		/*! \brief Numero di parole distinte
        */
		unsigned int W;
		
		/*! \brief Numero di parole totali
        */
		unsigned int N;
		
		/*! \brief Numero di documenti
        */
		unsigned int D;
		
		/*! \brief Contiene le etichette assegnate alle parole ad ogni iterazione
        */
		vector<vector<unsigned int>> Labels;    
 
		/*! \brief Matrice necessaria per individuare il least square clustering
		 *  L'elemento in posizione (i,j) è una stima Monte Carlo della probabilità che la parola i sia nello stesso topic della parola j
		 */
		deque<double> Pairwise_probabilities;
		
		/*! \brief Struttura in cui memorizzare i nomi dei documenti.
		 *  Il primo doc è quello con id 1 nel dataset, cioè ho la corrispondenza doc id- nome doc
		*/
		vector<std::string> docs;

	public:

		/*! \brief Costruttore di default
        */
		CategoricalPosteriorAnalysis() = default;
		
		/*! \brief Distruttore di default
        */
		~CategoricalPosteriorAnalysis()=default;

		/*! \brief Imposta la catena AllK
		 *  \return Lunghezza della catena
        */
		unsigned int SetAllK();
				
		/*! \brief Imposta la catena Alpha
		 *  \return Lunghezza della catena
        */
		unsigned int SetAllAlpha();
		
		/*! \brief Imposta la catena Gamma
		 *  \return Lunghezza della catena
        */
		unsigned int SetAllGamma();
		
		/*! \brief Imposta il vocabolario
        */
		//void SetVocabulary();
		
		/*! \brief Chiama lo script R che visualizza le traiettorie della distribuzione dei topic nel corpus
		 *  \param R - Istanza di R
		 *  \param BestClustering - iterazione a cui è stato individuato il least square clustering
        */
		//void VisualizeBeta(RInside& R, unsigned long BestClustering);
		
		/*! \brief Chiama lo script R per l'analisi della catena dei K
		 *  \param R - Istanza di R
		 *  \param Burnin - numero di valori iniziali della catena da scartare
		 *  \param Thinning - tiene un valore della catena ogni thinning valori
		*/
		void KPosteriorAnalysis(RInside& R, const unsigned int Burnin, const unsigned int Thinning);
		
		/*! \brief Chiama lo script R per l'analisi delle catene Alpha e Gamma
		 *  \param R - Istanza di R
		 *  \param AlphaBurnin - numero di valori iniziali della catena Alpha da scartare
		 *  \param AlphaThinning - tiene un valore della catena Alpha ogni AlphaThinning valori
		 *  \param GammaBurnin - numero di valori iniziali della catena Gamma da scartare
		 *  \param GammaThinning - tiene un valore della catena Gamma ogni GammaThinning valori
		 *  \param  AlphaTry - yes se si vuole ripetere l'analisi della catena Alpha, no altrimenti
		 *  \param  GammaTry - yes se si vuole ripetere l'analisi della catena Gamma, no altrimenti
		 */
		void AGPosteriorAnalysis(RInside& R, const unsigned int AlphaBurnin,const unsigned int AlphaThinning,const unsigned int GammaBurnin,
                                 const unsigned int GammaThinning, const char AlphaTry,const char GammaTry);
		
		/*! \brief Imposta la working directory di R
		 * \param _wd - working directory di R
		*/
		void Setwd(const std::string& _wd);
		
		/*! \brief Imposta il numero di parole distinte
		 * \param _W - numero di parole distinte
		*/
		void SetW(const unsigned int _W);
		
		/*! \brief Imposta il numero di documenti
		*   \param _D - numero di documenti
		*/
		void SetD(const unsigned int _D);
		
		/*! \brief Imposta il numero totale di parole 
		*   \param _N - numero totale di parole 
		*/
		void SetN(const unsigned int _N);
		
		/*! \brief Imposta la ricerca dei cluster, acquisendo da file i \f$ \theta \f$ e i  \f$ \beta \f$ delle ultime iterazioni
		 *  \param MaxIt - numero di iterazioni dell'algoritmo
		*/
		
		void SetTheta(const unsigned long MaxIt);
		
		/*! \brief Legge da file e imposta i \f$ \pi \f$
		 *  \param MaxIt - numero di iterazioni dell'algoritmo
		*/
		
		void SetPi(const unsigned long MaxIt);
		
		/*! \brief Acquisisce i nomi dei documenti
		*/
		void SetDocs();

		/*! \brief Acquisisce le etichette assegnate alle parole nelle ultime iterazioni
		 *  \param NrClusterings - numero di iterazioni da monitorare
		*/
		void LoadLabels(const unsigned long NrClusterings);  
		
		/*! \brief Individua il miglior clustering secondo il criterio dei minimi quadrati; si veda equazione (4.3) in Relazione_Parisi_Perego.pdf
		*   \return Iterazione a cui è stato individuato il miglior clustering
		*/
		unsigned long LeastSquareClustering();

		/*! \brief Calcola e stampa a video l'indice LPML; si vedano equazioni (4.1)- (4.2) in Relazione_Parisi_Perego.pdf
		*   \param MaxIt - numero di iterazioni da considerare nel calcolo di LPML
		*/
		//void LPML(unsigned long MaxIt);
		
		/*! \brief Effettua il riconoscimento dei topic
		*   \param R - Istanza di R
		*/
		//void TrackingClusters(RInside& R);
		
		/*! \brief Associa i documenti ai topic
		*/
		
		void AssociatingDocs();
		

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
template<unsigned int DIM> unsigned int CategoricalPosteriorAnalysis<DIM>::SetAllK(){

    std::ifstream file("./cpp_results/AllK.txt");

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

	return AllK.size();

}

//Imposta la catena Alpha
template<unsigned int DIM> unsigned int CategoricalPosteriorAnalysis<DIM>::SetAllAlpha(){

    std::ifstream file ("./cpp_results/AllAlpha.txt");

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

	return Alpha.size();

}

//Imposta la catena AllK
template<unsigned int DIM> unsigned int CategoricalPosteriorAnalysis<DIM>::SetAllGamma(){

    std::ifstream file ("./cpp_results/AllGamma.txt");

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

	return Gamma.size();

}

/*
//Imposta il vocabolario
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetVocabulary(){

    std::ifstream file("Vocabulary.txt");

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
*/
/*
//Chiama lo script R che visualizza le traiettorie della distribuzione dei cluster nel corpus
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::VisualizeBeta(RInside& R, unsigned long BestClustering){

    try {
		//Devo andare nella cartella base
        R.parseEvalQ(wd);

		//unordered_map<unsigned int,vector<tuple<TopicId,Check,double,vector<double>>>> Theta;

		unsigned int zeros = 0;
		unsigned int maxId = 0;
		unsigned int tmp_maxId = 0;

		for(size_t it = 0; it< Theta.size() ; ++it){

			for(size_t k=0; k< Theta[it].size(); ++k){

				if(std::get<0> (Theta[it][k]) == 0)
					++ zeros;

				tmp_maxId = std::get<0> (Theta[it][k]);
				maxId = std::max(maxId,tmp_maxId);

			}

		}


		AllBeta.reserve(Theta.size());
		vector<double> tmp_beta( maxId+zeros ,0.0 );
		size_t first_pos = 0; // la prima posizione non vuota in cui piazzare il peso del cluster con etichetta 0

		for(size_t it = 0; it< Theta.size() ; ++it){
			AllBeta.push_back(tmp_beta);
			for(auto i : Theta[it]){
				if(std::get<0>(i) > 0)
					AllBeta[it][std::get<0>(i) + zeros - 1] = std::get<2>(i); //beta
				else{
					AllBeta[it][first_pos] = std::get<2>(i); //beta
					++ first_pos;
				}
			}
		}
		
		//calcolo la media: ogni vettore in AllBeta ha la stessa dimensione. 
		//Medio componente per componente, questa è La traiettoria media 
		
		vector<double> BetaMean;
		double temp_beta_mean = 0.0;
		BetaMean.reserve(maxId+zeros);
		
		
		for(size_t m = 0; m< maxId+zeros ; ++m){
			for(size_t it = 0; it< AllBeta.size() ; ++it){
				temp_beta_mean += AllBeta[it][m];
			}
			
			temp_beta_mean = temp_beta_mean/static_cast<double>(AllBeta.size());
			
			BetaMean.push_back(temp_beta_mean);
		}


		Rcpp::List RBeta(AllBeta.size());

		double temp_weight = 0.0;

		for(size_t k= 0; k<AllBeta.size(); ++k){

			for(size_t i = 0; i< AllBeta[k].size(); ++i)
				temp_weight += AllBeta[k][i];

			AllBeta[k].push_back(1 - temp_weight);
			RBeta[k]= AllBeta[k];
			temp_weight = 0.0;
		}
		
		Rcpp::NumericVector RBetaMean(BetaMean.size());
		
		for(size_t m = 0; m< BetaMean.size(); ++m){
			RBetaMean[m]= BetaMean[m];
		}
		
		// Ora posso mandare il vettore nel workspace di R

		if(BestClustering < std::numeric_limits<unsigned long>::max())
			R["BestClustering"] = BestClustering + 1;   // R comincia a contare da 1
		else
			R["BestClustering"] = -1;
		

		vector<unsigned int> labels;
		vector<unsigned int> at;
		labels.reserve(zeros+maxId+1);
		at.reserve(zeros+maxId+1);


		for(size_t z=0;z<zeros;++z)
			labels.push_back(0);

		for(size_t l=0 ; l<= maxId; ++l)
			labels.push_back(l+1);

		for(size_t t=0; t<labels.size() ; ++t)
			at.push_back(t);

		Rcpp::NumericVector Rlabels(labels.size());
		Rcpp::NumericVector Rat(at.size());

		for(size_t l = 0; l< labels.size() ; ++l)
			Rlabels[l] = labels[l];

		for(size_t a = 0; a< at.size() ; ++a)
			Rat[a] = at[a];


        R["Beta"] = RBeta;
		R["BetaMean"] = RBetaMean;
        R["Dim"] = AllBeta.size();
		R["TopicLabels"] = Rlabels;
		R["at"] = Rat;

        // Chiamo lo script
        std::string src = "source(\"../../src/Rscript/Beta.R\")";
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
*/

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
	return;
}


//Acquisisce le etichette assegnate ai dati nelle ultime iterazioni
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::LoadLabels(const unsigned long NrClusterings){

	Labels.resize(NrClusterings);

	std::ifstream file("./cpp_results/Labels.bin",std::ios::binary);

	if(file.fail()){
		std::cerr<<"Error in LoadLabels: cannot open Labels.bin"<<std::endl;
		exit(1);
	}

	for(size_t it = 0; it < NrClusterings; ++it){

		unsigned int* tmp = new unsigned int[N];

		file.read(reinterpret_cast<char*>(tmp), N*sizeof(unsigned int));

		Labels[it].assign(tmp,tmp + N);

		delete [] tmp;
	}


    file.close();

}


//Individua il miglior clustering secondo il criterio dei minimi quadrati
template<unsigned int DIM> unsigned long CategoricalPosteriorAnalysis<DIM>:: LeastSquareClustering(){

	unsigned long long dim = (N*(N+1))/2;

	Pairwise_probabilities.resize(dim);

	//inizializzo il vettore con zeri
	for(size_t i = 0; i< dim; ++i)
		Pairwise_probabilities[i] = 0.0;

	unsigned long long indices;

	//metto 1 sulla diagonale : il dato i-esimo ha probabilità 1 di finire nello stesso cluster di se stesso
	for(size_t i = 0;i<N;++i){
		indices = ((i+1)*i)/2 +i;
		Pairwise_probabilities[indices] = 1.0;
	}

	// calcolo le frequenze assolute degli eventi (i ha la stessa etichetta di j) per ogni i e j

	for(size_t it = 0; it < Labels.size(); ++it){
		for(size_t i = 0; i<N; ++i){
			for(size_t j = 0; j<i; ++j){
				indices = ((i+1)*i)/2 +j;
				if(Labels[it][i] == Labels[it][j]) Pairwise_probabilities[indices]++;
			}
		}
	}

	// calcolo le frequenze relative -> stime delle pairwise probability

	for(size_t i = 0; i<N; ++i){
		for(size_t j = 0; j<i; ++j){
			indices = ((i+1)*i)/2 +j;
			Pairwise_probabilities[indices] /= static_cast<double>(Labels.size());
		}
	}


	double minDB = std::numeric_limits<double>::max();     // il massimo double possibile
	unsigned long minIndex = 0;
	double db;
	double deviance;

	for(size_t it = 0; it<Labels.size(); ++it){
		db = 0.0;
		for(size_t i = 0; i<N; ++i){
			for( size_t j=0; j<i; ++j){
				indices = ((i+1)*i)/2 +j;
				deviance = (Labels[it][i] == Labels[it][j] ? 1 : 0 ) - Pairwise_probabilities[indices];
				db += deviance*deviance;
			}
		}

		if( db < minDB ) {              //devo trovare il clustering che minimizza db
			minDB = db;                     // questo è il minimo
			minIndex = it;                   // questo è il clustering che realizza il minimo
		}
	}

	unsigned int K = *(std::max_element(Labels[minIndex].cbegin(), Labels[minIndex].cend()))+1;

	vector<unsigned int> BestClustering(K,0);

	for(std::vector<unsigned int>::const_iterator it = Labels[minIndex].cbegin(); it != Labels[minIndex].cend(); ++it)
		++BestClustering[*it];

	std::cout<<"Best clustering found at iteration "<<minIndex<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Best clustering has "<<K<<" clusters"<<std::endl;
	std::cout<<std::endl;
	for(unsigned int k = 0; k<K;++k)
		std::cout<<"Cluster "<<k<<" has "<<BestClustering[k]<<" elements"<<std::endl;
	std::cout<<std::endl;
	
	// ho individuato l'iterazione in cui c'è la partizione ottima (minIndex), ora recupero beta, pi, theta di quella iterazione
	// e li memorizzo negli attributi appositi
	
	Pi_opt = AllPi[minIndex];
	
	vector<tuple<GroupId,Check,double,vector<double>>> q_opt = Theta[minIndex];
	
	for(auto iter = q_opt.cbegin(); iter!= q_opt.cend(); ++iter){
		
		Theta_opt.push_back(std::get<3>(*iter));
		Beta_opt.push_back(std::get<2>(*iter));
		
	}

	for(auto i : Pi_opt){
		for(auto j : i){
			std::cout<<j<<" ";
		}
		std::cout<<std::endl;
	}

	std::cout<<std::endl;

	for(auto i : Beta_opt)
		std::cout<<i<<" ";
	std::cout<<std::endl;
    std::cout<<std::endl;

	for(auto i: Theta_opt){
		for(auto j: i)
			std::cout<<j<<" ";
	std::cout<<std::endl;


	}

	
	return minIndex;

}


//Imposta la ricerca dei cluster, acquisendo da file i \theta e i \beta delle ultime iterazioni
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetPi(const unsigned long MaxIt){
	
	std::ifstream LastPi("./cpp_results/LastPi.bin",std::ios::binary);

    if(LastPi.fail()){
		std::cerr<<"Cannot open file: LastPi"<<std::endl;
		exit(1);
    }
	
	unsigned int K;
	
	for(size_t it=0; it< MaxIt; ++it){
		
		K = AllK[it+1];    // perchè ho scartato il primo K
		
		for(size_t d=0; d<D; ++d){
			
			double *tmp_pi = new double [K];
			LastPi.read(reinterpret_cast<char*>(tmp_pi), (K)*sizeof(double));
            std::vector<double> TempPi(tmp_pi, tmp_pi + (K));
            delete [] tmp_pi;
			
			AllPi[it].push_back(TempPi);
		}
				
	}

	/*
	std::ofstream file2("PiPost.txt");

	for(size_t it = 0; it< MaxIt; ++it){
		for(auto iter = AllPi[it].cbegin(); iter != AllPi[it].cend(); ++iter){
            for(auto i : *iter)
				file2<<i<<" ";
			file2<<std::endl;
		}
		file2<<"###########"<<std::endl;
	}
	*/

	
	
}



//Acquisisce i nomi dei documenti
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetDocs(){
	
	std::ifstream file("Docs.txt");
	
	if(file.fail()){
		std::cerr<<"Cannot open file: Docs.txt"<<std::endl;
		exit(1);
    }
	
	std::string ss;
	
	while(std::getline(file,ss))
		docs.push_back(ss);
	
	if(docs.size() != D){
		std::cerr<<" read nr of documents different from D"<<std::endl;
		exit(1);
    }

	/*for(auto i : docs)
		std::cout<<i<<std::endl;
	*/
}


//Imposta la ricerca dei cluster, acquisendo da file i \theta e i \beta delle ultime iterazioni
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::SetTheta(const unsigned long MaxIt){

    std::ifstream LastTheta("./cpp_results/LastTheta.bin",std::ios::binary);

    if(LastTheta.fail()){
		std::cerr<<"Doesn't open file: LastTheta"<<std::endl;
		exit(1);
    }

    std::ifstream LastBeta("./cpp_results/LastBeta.bin",std::ios::binary);

    if(LastBeta.fail()){
		std::cerr<<"Doesn't open file: LastBeta"<<std::endl;
		exit(1);
    }

    bool check = true;
    vector<tuple<GroupId,Check,double,vector<double>>> TempTheta;
    unsigned int iter = 0;


    //while(check){

      //  unsigned int ITER_Theta;
       // LastTheta.read(reinterpret_cast<char*>(&ITER_Theta), sizeof(unsigned int));
       //unsigned int ITER_Beta;
       // LastBeta.read(reinterpret_cast<char*>(&ITER_Beta), sizeof(unsigned int));

	    unsigned int K;
	   
        for(size_t it = 0; it<MaxIt; it++){

            K = AllK[it+1] ;
            //LastTheta.read(reinterpret_cast<char*>(&K), sizeof(unsigned int));
            //unsigned int d;
            //LastBeta.read(reinterpret_cast<char*>(&d), sizeof(unsigned int));
            double *tmp_beta = new double [K];
            LastBeta.read(reinterpret_cast<char*>(tmp_beta), K*sizeof(double));
            std::vector<double> TempBeta(tmp_beta, tmp_beta + K);
            delete [] tmp_beta;

            for(size_t k=0; k<K; ++k){
                double *tmp_theta=new double[W];
                LastTheta.read(reinterpret_cast<char*>(tmp_theta),W*sizeof(double));
                std::vector<double> fun(tmp_theta, tmp_theta + W);
                delete [] tmp_theta;
                TempTheta.push_back(std::make_tuple(0,0,TempBeta[k],fun));
            }

            Theta.insert({iter,TempTheta});
            iter ++;
            TempTheta.clear();
        } //chiudo for

		//check=false;
   // } // chiudo while

    LastBeta.close();
    LastTheta.close();

	/*
	std::ofstream file1("BetaPost.txt");

	for(size_t it = 0; it< MaxIt; ++it){
		for(auto iter = Theta[it].cbegin(); iter != Theta[it].cend(); ++iter)
            file1<<std::get<2>(*iter)<<" ";
		file1<<std::endl;
		
		file1<<"###########"<<std::endl;
	}

	std::ofstream file2("ThetaPost.txt");

	for(size_t it = 0; it< MaxIt; ++it){
		for(auto iter = Theta[it].cbegin(); iter != Theta[it].cend(); ++iter){
            for(auto iter_in = std::get<3>(*iter).cbegin(); iter_in != std::get<3>(*iter).cend(); ++iter_in)
				file2<<*iter_in<<" ";
			file2<<std::endl;
		}
		file2<<"###########"<<std::endl;
	}
	*/

}

/*
//Associa i documenti ai topic
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::AssociatingDocs(){
	
	unsigned int K_opt = Theta_opt.size();
	
	//mi interessa classificare i topics in base ai beta? qui per ora no, forse mi serve per mostrare i wordcloud
	
	//devo avere la corrispondenza Id doc - nome doc
	
	// per ogni j, cerco in \pi_j il peso massimo \pi_jk e associo al doc j il topic k-1
	// nel vettore corrisp, lungo K, l'elemento k è un vettore di doc associati al topic k,
	// se il vettore è vuoto vuol dire che a quel topic non ho associato nessun doc (controllo!)
	// in questo modo ogni doc è stato associato ad uno e un solo topic
	// alla fine per ogni topic k stampo i nomi dei doc che gli sono stati associati
	
	vector<vector<std::string>> topic_to_docs(K_opt);
	
	for(size_t d = 0; d< Pi_opt.size(); ++d){
		auto max_pi = std::max_element(Pi_opt[d].cbegin(),Pi_opt[d].cend());
		unsigned int position = std::distance(Pi_opt[d].cbegin(),max_pi);
		if(position != K_opt) //vuole dire che il peso più grande è quello del cluster vuoto, non dovrebbe succedere
			topic_to_docs[position].push_back(docs[d]);
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
		}
		
	}
		
}
*/


//Associa i documenti ai topic_ Versione 2
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::AssociatingDocs(){
	
	// per ogni topic k, controllo per ogni doc j se \pi_jk è maggiore di una certa soglia
	// per il momento fisso la soglia, poi verrà decisa dall'utente e potrebbe variare per ogni topic (interazione utente)
	// se la condizione della soglia è vera, inserisco l'id del documento nel vettore topics_to_docs, poi è tutto come prima
	// con questo metodo non ho la certezza di associare ogni doc a qualche topic, perciò metto un flag che diventa vero se associo un doc ad almeno un topic
	
	unsigned int K_opt = Theta_opt.size();	

	double threshold = 0.5;//1.0/K_opt;
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
	
	std::ofstream file("DocAnalysis.txt");
	
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
}

/*

//Calcola e stampa a video l'indice LPML
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::LPML(const unsigned long MaxIt){

    std::ifstream lpml ("./../cpp_results/CPO.bin");
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
*/
/*
//Effettua il riconoscimento dei cluster
template<unsigned int DIM> void CategoricalPosteriorAnalysis<DIM>::TrackingClusters(RInside& R){

    size_t iter = Theta.size();
    std::vector<double> theta_itk;
    std::vector<size_t> idx(W,0);
    for(size_t i=0; i<W; i++)
        idx[i] = i;

    unsigned int MinCheck = 0;
    unsigned int tmp_MinCheck;
    double sum = 0.0;
    //decidiamo il numero di parole minimo per far si che i cluster siano rappresentati almeno per il loro 70%
    for(size_t It=0; It<iter; It++){
        for(auto ClustersIter = Theta[It].begin(); ClustersIter != Theta[It].end(); ClustersIter++){
            theta_itk = std::get<3>(*ClustersIter);
            //rioridino gli indici di theta
            std::sort(idx.begin(),idx.end(),[&theta_itk](size_t i1, size_t i2) {return theta_itk[i1] > theta_itk[i2];});
            //somma dei pesi fin quando non raggiungo 0.7
            sum = 0.0;
            tmp_MinCheck = 0;
            for(size_t i=0; i<W;i++){
                sum += theta_itk[idx[i]];
                tmp_MinCheck++;
                if(sum > 0.9) break;
            }
            if(tmp_MinCheck > MinCheck)
                MinCheck = tmp_MinCheck;
        }
    }


    for(size_t It=0; It<iter; It++){
        for(auto ClustersIter = Theta[It].begin(); ClustersIter != Theta[It].end(); ClustersIter++){
            theta_itk = std::get<3>(*ClustersIter);
            //rioridino gli indici di theta
            std::sort(idx.begin(),idx.end(),[&theta_itk](size_t i1, size_t i2) {return theta_itk[i1] < theta_itk[i2];});
            //azzeramento pesi che non ci interessano
            for(size_t i=0; i<W-MinCheck;i++)
                theta_itk[idx[i]] = 0.0;
            //riassegnamo al Theta
            std::get<3>(*ClustersIter) = theta_itk;

        }
    }
    std::cout<<"**********************************************"<<std::endl;
    std::cout<<"             Research of Clusters             "<<std::endl;
    std::cout<<"**********************************************"<<std::endl;
    std::cout<<std::endl;
	std::cout<<"Sorry, it will take a while..."<<std::endl;
	std::cout<<std::endl;
    //inizio a fare la ricerca
    using namespace CoeffSimilitudine;

    double coef;
    double tmp_coef;
    std::vector<double> theta_conf;
    std::vector<double> tmp_theta;
    GroupId topic = 1;
    Check nwords, tmp_nwords;
    double tmp_sum = 0.0;

    //blocco di appartenenza,posizione nel blocco, conteggio, vettore per il controllo
    std::vector<std::tuple<size_t,size_t,Check>> TempCluster;
    std::vector<std::tuple<GroupId,Check,double,std::vector<double>>> RifBlock;
    std::vector<std::tuple<GroupId,Check,double,std::vector<double>>> ConfBlock;

    size_t block_pos,vect_pos,tmp_pos;

    for(size_t It =0; It<iter - 1; It++){

     //   std::cout<<"---------------BLOCCO RIFERIMENTO  "<<It<<"----------------------"<<std::endl;
        RifBlock = Theta[It];

        for(size_t It_RifBlock = 0; It_RifBlock<RifBlock.size(); It_RifBlock++){

             if(std::get<0>(RifBlock[It_RifBlock]) == 0){
                 theta_itk = std::get<3>(RifBlock[It_RifBlock]);
                 TempCluster.clear();
                 TempCluster.push_back(std::make_tuple(It,It_RifBlock,0));

                 for(size_t It_conf=It+1; It_conf<iter; It_conf++){
                  //  std::cout<<"---------------BLOCCO CONFRONTO  "<<It_conf<<"----------------------"<<std::endl;
                    ConfBlock = Theta[It_conf];
                    coef = 0.0;
                    nwords = 0;
                    tmp_pos=0;
                    for(size_t It_ConfBlock=0; It_ConfBlock<ConfBlock.size(); It_ConfBlock++){
                        tmp_coef = Coeff(theta_itk,std::get<3>(ConfBlock[It_ConfBlock]));

                         if(coef < tmp_coef){
                            tmp_theta = std::get<3>(ConfBlock[It_ConfBlock]);
                            tmp_nwords =0;
                            tmp_sum = 0.0;

                            for(size_t w=0; w<W; w++){
                                if(theta_itk[w]!=0.0 && tmp_theta[w]!=0.0){
                                    tmp_nwords++;
                                    tmp_sum += tmp_theta[w];
                                }
                            }
                            if(tmp_sum > 0.69 && std::get<1>(ConfBlock[It_ConfBlock])<tmp_nwords && nwords<tmp_nwords){
                                nwords = tmp_nwords;
                                coef = tmp_coef;
                                tmp_pos = It_ConfBlock;
                            }

                         }

                    }

                    TempCluster.push_back(std::make_tuple(It_conf,tmp_pos,nwords));
                 }
                //std::cout<<"Dim: "<<TempCluster.size()<<std::endl;
                //assegnamo i cluster e i check
                if(TempCluster.size()>1){

                    for(size_t It_TC =0; It_TC<TempCluster.size() - 1; It_TC ++){
                        ConfBlock = Theta[std::get<0>(TempCluster[It_TC])];
                        vect_pos = std::get<1>(TempCluster[It_TC]);
                        theta_conf = std::get<3>(ConfBlock[vect_pos]);

                        for(size_t Itconf_TC = It_TC +1; Itconf_TC<TempCluster.size(); Itconf_TC++){
                            ConfBlock = Theta[std::get<0>(TempCluster[Itconf_TC])];
                            vect_pos = std::get<1>(TempCluster[Itconf_TC]);
                            tmp_theta = std::get<3>(ConfBlock[vect_pos]);

                            tmp_nwords = 0;
                            for(size_t w=0; w<W; w++){
                                if(theta_conf[w]!=0.0 && tmp_theta[w]!=0.0)
                                    tmp_nwords++;
                            }
                            block_pos = std::get<0>(TempCluster[It_TC]);
                            vect_pos = std::get<1>(TempCluster[It_TC]);
                            if(std::get<1>(Theta[block_pos][vect_pos]) < tmp_nwords)
                                std::get<1>(Theta[block_pos][vect_pos]) = tmp_nwords;

                            block_pos = std::get<0>(TempCluster[Itconf_TC]);
                            vect_pos = std::get<1>(TempCluster[Itconf_TC]);
                            if(std::get<1>(Theta[block_pos][vect_pos]) < tmp_nwords)
                                std::get<1>(Theta[block_pos][vect_pos]) = tmp_nwords;

                        }

                    }
                    for(size_t It_TC =0; It_TC<TempCluster.size(); It_TC ++){
                        block_pos = std::get<0>(TempCluster[It_TC]);
                        vect_pos = std::get<1>(TempCluster[It_TC]);

                        std::get<0>(Theta[block_pos][vect_pos]) = topic;
                    }
                    topic++;

                }

             }
        }
    }

    topic--;
    //individuiamo la calssifica dei topic più frequenti
    vector<unsigned int> Frequency(topic,0);
    unsigned int n;
    char Try = 'y';
    GroupId k;
    for(size_t It =0; It<iter; It++){
        RifBlock = Theta[It];
        for(size_t It_RifBlock=0; It_RifBlock<RifBlock.size(); It_RifBlock++){
            k = std::get<0>(RifBlock[It_RifBlock]);
            if(k!=0)
                Frequency[k-1]++;
        }
    }
    //oridno in modo decrescente
    idx.clear();
    for(size_t i =0; i<topic; i++)
        idx.push_back(i);


    std::sort(idx.begin(),idx.end(),[&Frequency](size_t i1, size_t i2) {return Frequency[i1] > Frequency[i2];});
    std::cout<<"Topic found: "<<topic<<" labels(1,...,"<<topic<<")"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"The Topic ordered from the most frequent to the less frequent: ";
    for(auto i: idx)
        std::cout<<i+1<<" ";
    std::cout<<std::endl;

    vector<vector<double>> AllKBeta;
    vector<double> tmp_Beta;
    bool check_beta;


    while(Try == 'y'){
        std::cout<<std::endl;
        std::cout<<"How many topics do you want visualize?  n = ";
        std::cin>>n;

        while(n>topic){
            std::cout<<"insert again a number between 1 - "<<topic<<": ";
            std::cin>>n;
        }

        AllKBeta.clear();
        check_beta = false;

        for(size_t It_idx=0; It_idx<static_cast<size_t>(n); It_idx++){
           //individuo il blocco
            tmp_Beta.clear();
            std::cout<<"Searching topic "<<idx[It_idx] +1<<std::endl;
            for(size_t It =0; It<iter; It++){

                RifBlock = Theta[It];

                for(size_t It_RifBlock = 0; It_RifBlock<RifBlock.size(); It_RifBlock++){

                    if(std::get<0>(RifBlock[It_RifBlock]) == idx[It_idx] + 1){
                        //std::cout<<"Topic "<<std::get<0>(RifBlock[It_RifBlock])<<std::endl;
                        tmp_Beta.push_back(std::get<2>(RifBlock[It_RifBlock]));
                        check_beta = true;
                        break;
                    }
                    else check_beta = false;
                }
                if(!check_beta)
                tmp_Beta.push_back(0.0);
            }
            AllKBeta.push_back(tmp_Beta);
        }
        try {
		//Devo andare nella cartella base
        R.parseEvalQ(wd);

        // Trasformo Beta dalla codifica in c++ a NumericVector, importabile in R
        //const unsigned int K = AllBeta.size();
		Rcpp::List RBeta(AllKBeta.size());
		Rcpp::NumericVector RTopics(n);

		for(size_t k= 0; k<AllKBeta.size(); ++k)
			RBeta[k]= AllKBeta[k];
        for(size_t It_idx = 0; It_idx<n; It_idx++)
            RTopics[It_idx] = idx[It_idx] + 1;
		// Ora posso mandare il vettore nel workspace di R

        R["Beta"] = RBeta;
        R["Dim"] = AllKBeta.size();
        R["Labels"] = RTopics;

        // Chiamo lo script
        std::string src = "source(\"../../src/Rscript/TrackingClusters.R\")";
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

        std::cout<<"Would you like change the number of topic to visualize? press y [yes] or n [no] ";
        std::cin>>Try;

        while(Try !='y' && Try != 'n'){
            std::cout<<"press y [yes] or n [no] "<<std::endl;
            std::cin>>Try;
        }

    }

    //WordCloud

    Try='y';
    unsigned int w;
    double tmp_MaxBeta;
    vector<vector<double>> weights;
    vector<vector<unsigned int>> WCW;
    vector<unsigned int> tmp_WCW;
    vector<vector<DataId>> Wordclouds;
    vector<DataId> tmp_Wordclouds;
    vector<size_t> idx_word;
    GroupId tmp_k;

    for(size_t i=0; i<W; i++)
        idx_word.push_back(i);


    std::cout<<std::endl;
    std::cout<<"**********************************************"<<std::endl;
    std::cout<<"             Creating WordClouds               "<<std::endl;
    std::cout<<"**********************************************"<<std::endl;
    std::cout<<std::endl;
    while(Try == 'y'){
        std::cout<<std::endl;
        std::cout<<"How many topics do you want to visualize?  n = ";
        std::cin>>n;

        while(n>topic){
            std::cout<<"insert again a number between 1 - "<<topic<<": ";
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
        //cerco i cluster scelti
        weights.clear();
        for(size_t It_idx=0; It_idx<static_cast<size_t>(n); It_idx++){
           //individuo il blocco
            tmp_k =idx[It_idx] +1;
            std::cout<<"Searching topic "<<tmp_k<<" ..."<<std::endl;
            tmp_MaxBeta = 0.0;
            for(size_t It =0; It<iter; It++){

                RifBlock = Theta[It];
                for(size_t It_RifBlock = 0; It_RifBlock<RifBlock.size(); It_RifBlock++){
                    if(std::get<0>(RifBlock[It_RifBlock])==tmp_k && std::get<2>(RifBlock[It_RifBlock])>tmp_MaxBeta){
                        tmp_MaxBeta = std::get<2>(RifBlock[It_RifBlock]);
                        tmp_theta = std::get<3>(RifBlock[It_RifBlock]);
                    }
                }
            }
            weights.push_back(tmp_theta);
        }
        std::cout<<std::endl;
        //prendo le parole che mi interessano
        WCW.clear();
        Wordclouds.clear();
        for(size_t it_weight=0; it_weight<weights.size(); it_weight++){
            tmp_theta = weights[it_weight];

            std::sort(idx_word.begin(),idx_word.end(),[&tmp_theta](size_t i1, size_t i2) {return tmp_theta[i1] > tmp_theta[i2];});
            tmp_WCW.clear();
            tmp_Wordclouds.clear();
            for(size_t i=0; i<w; i++){
                tmp_WCW.push_back(static_cast<unsigned int>(tmp_theta[idx_word[i]]*N));
                tmp_Wordclouds.push_back(idx_word[i] + 1);
            }
            WCW.push_back(tmp_WCW);
            Wordclouds.push_back(tmp_Wordclouds);

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
		Rcpp::List RWeights(WCW.size());
		Rcpp::List RWords(Wordclouds.size());
		Rcpp::CharacterVector RVocab(W);
		Rcpp::NumericVector RTopics(n);

		for(size_t k= 0; k<WCW.size(); ++k)
			RWeights[k]= WCW[k];
        for(size_t k= 0; k<Wordclouds.size(); ++k)
			RWords[k]= Wordclouds[k];

        for(size_t It_idx = 0; It_idx<n; It_idx++)
            RTopics[It_idx] = idx[It_idx] + 1;

        for(unsigned int i=1; i<=W; i++)
            RVocab[i-1] = Vocabulary[i];

		// Ora posso mandare il vettore nel workspace di R

        R["Weights"] = RWeights;
        R["Dim"] = WCW.size();
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
*/
#endif

