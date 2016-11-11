#ifndef __HDP_MCMC__HPP__
#define __HDP_MCMC__HPP__

#include "Model.hpp"
#include "Document.hpp"
#include "Struct.hpp"
#include "omprng.hpp"
#include "Functions.hpp"

#include <random>
#include <tuple>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <omp.h>


using std::vector;
using std::log;

/*! \mainpage Hierarchical Dirichlet Process Mixture Model parallelizzabile
 *
 *  \section intro Introduzione
 *  In questo codice e' implementato l'algoritmo proposti nel paper:
 *  J. Chang, J. W. Fisher III Parallel Sampling of HDPs using sub-clusters splits,
 *  NIPS, 2014.
 *
 *  Nel paper si fa riferimento solo al problema del topic-modeling, ma l'algoritmo è estendibile a problemi di altro tipo, vedere il capitolo 7
 *  di Relazione_Parisi_Perego.pdf
 *  In questo algoritmo si alternano passi di Gibbs sampler a passi di Metropolis-Hastings.
 *  Nei passi di Gibbs sampler si considerano soltanto i cluster non vuoti e si aggiornano le quantità di interesse \f$ (\beta,\pi,\theta) \f$
 *  usando le full-conditional, equazioni dalla 3.3 alla 3.11;
 *  queste ultime non sono altro che alcuni step dell'algoritmo. Durante questi passi il numero di cluster rimane invariato.
 *  Con i passi di Metropolis-Hastings, invece, si propone l'unione di due cluster (mosse di merge) oppure la divisione di un cluster (mosse di split).
 *  Per ogni topic \f$ k \f$ si individuano due sub-topic, \f$ kl \f$ e \f$ kr \f$ che rispettivamente corrispondono al sub-topic sinistro
 *  e destro; i nuovi topic si propongono sulla base di tali sub-topic.
 *
 *  Per quanto riguarda l'inferenza sui parametri latenti e' stato trattato solo il caso di modelli
 *  coniugati.
 */

/*! \file HDP_MCMC.hpp
 *
 *  Classe di HDP_MCMC per l'esecuzione dell'algoritmo
 * \date  Febbraio 2016
 */

/*! \brief HDP_MCMC
 *
 * Implementazione algoritmo globale.
 * Gibbs Sampler per i clusters e i sub-cluster, Metropolis-Hastings per le mosse di Merge/Split locale e globale.
 * Sono implementate le equazioni di campionamento dei pesi globali e dei tavoli.
 * Si occupa di aggiornare i conteggi dei cluster in Model, controllando la situazione dei gruppi.
 * Per dettagli sul funzionamento dell'algoritmo globale consultare documento
 * Relazione_Parisi_Perego.pdf capitolo 3.
 * Questa classe si appoggia sulle classi Cluster.hpp per la gestione dei cluster e la definizione di verosimiglianza,
 * Model.hpp  per l'eventuale campionamento dalle distribuzioni di interesse, l'inferenza sui parametri
 * latenti e la definizione di prior e invece si appoggia a Document.hpp per la gestione dei dati e per il campionamento delle etichette
 * La particolare classe Model e' l'istanzazione del parametro template MODEL.
 * In Model.hpp e' implementato CategoricaModel (prior Dirichlet).
 * \date Febbrario 2016
 */
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT, unsigned int DIM=1>
class HDP_MCMC final
{
    public:

        using THETA = typename MODEL<DIM>::THETA;
        using POINT = typename MODEL<DIM>::POINT;
        using HYP = typename MODEL<DIM>::HYP;
        using Corpus = vector<DOCUMENT<DIM>>;
        using ClusterID = unsigned int;

	private:

        /*! \brief Oggetto che gestisce i dati
         */
        Corpus corpus;

    /*! \brief Flags che fornisce le informazioni necessarie per eseguire l'algoritmo con le caratteristiche richieste
     *  dall'utente.
     *  Il primo elemento e' true se il dateset e' stato caricato.
     *  Il secondo elemento e' true se sono state inserite informazioni su Alpha
     *  Il terzo elemento e' true se si vuole eseguire l'algoritmo con le prior su Alpha
     *  Il terzo elemento è false se si vuole eseguire l'algoritmo con Alpha fisso
     *  Il quarto elemento è true se sono stati inserite informazioni su Gamma
     *  Il quinto elemento è true se si vuole eseguire l'algoritmo con le prior su Gamma
     *  Il quinto elemento è false se si vuole eseguire l'algortimo con Gamma fisso
     *  Il sesto elemento è true se si sono inserite informazioni su Lambda
     *  Il settimo elemento è true se si è inserito un numero iniziale di cluster
     *  L'ottavo elemento è true se si è scelto un seme
     *  Il non elemento è true se si vuole monitorare il modello con LPML.
     */
        std::tuple<bool,bool,bool,bool,bool,bool,bool,bool,bool> Flags;

        /*! \brief Numero totale di gruppi
         */
        unsigned int D;

        /*! \brief Numero totale di dati
         */
        unsigned int N;

        /*! \brief Dimensione dell'iperparametro della distribuzione del parametro latente
         */
        unsigned int W;

        /*! \brief Iterazione corrente
         */
        unsigned long It;

        /*! \brief Numeror massimo iterazioni (criterio stop)
         */
        unsigned long MaxIt;

        /*! \brief Numero di iterazioni per i sub-cluster (criterio di stop)
         */
        unsigned long MaxIt_SubCluster;

        /*! \brief  Parametro di concentrazione del processo di dirchlet che governa i clusters (all'iteriazione It)
         */
        double Gamma;

        /*! \brief Storico dei Gamma in caso dell'uso della prior.
         */
        vector<double>AllGamma;

        /*! \brief Parametro di forma della prior su Gamma
         */
        double GA;

        /*! \brief Parametro di rate della prior su Gamma
         */
        double GB;

        /*! \brief  Parametro di concentrazione del processo di dirchlet che governa i cluster nel gruppo (all'iteriazione It)
         */
        double Alpha;

        /*! \brief Strico degli Alpha nel caso dell'uso della prior
         */
        vector<double>AllAlpha;

        /*! \brief Parametro di forma della prior su Alpha
         */
        double AA;

        /*! \brief Parametro di rate della prior su Alpha
         */
        double AB;

        /*! \brief Numero corrente di cluster
         */
        unsigned int K;

        /*! \brief Storico dei K
         */
        vector<unsigned int> AllK;

        /*! \brief Stirling Number in scala logaritmica
         */
        vector<long double> LogStirlingNumbers;

        /*! \brief	Modello per la distribuzione iniziale.
		 *
		 * Tiene conto dell'opportuna verosimiglianza per
		 * avere un modello bayesiano coniugato per quanto riguarda l'inferenza sui parametri
		 * latenti.
		 */
        MODEL<DIM> Model;

        /*! \brief peso globale del cluster "vuoto" all'iterazione It.
         *  aggrega i pesi dei cluster (al momento vuoti), che potrebbero manifestarsi nelle iterazioni
         *  successive alla It-esima.
         */
        double Beta_empty;

        /*! \brief Matrice necessaria per il calcolo dell' Hasting Ratio
         * calcolata con le equazioni 8.12 e 8.13 di Relazione_Parisi_Perego.pdf
         */
        vector<double> logL;

        /*! \brief Numero totale di tavoli, da usare per campionare la Gamma, viene aggiornato dopo aver campionato
         *  dall'equazione 3.3 di Relazione_Parisi_Perego.pdf
         */
        unsigned int m;

        /*! \brief Oggetto che memorizza \f$ \sum_{It=1}^{MaxIt} f_{ij}(y_{ij}|\theta_{z_{ij}}^{(g)}) \f$
         */
        vector<double> CPO;

        /*! \brief Iterazione dal quale iniziare a calcolare il CPO
         */
        unsigned long burnin_CPO;

        /*! \brief Generatore di numeri casuali in parallelo
         */
        omprng Gen;

        /*! Numero di threads
         */
        unsigned int OMP_NUM_THREADS;


    //METODI PUBBLICI
	public:

        //COSTRUTTORE e DISTRUTTORE DI DEFAULT
        /*! \brief Costruttore di default della classe in cui tutti gli elementi vegono inzializzati con il loro costruttore di default.
         *  Ai valori scalari viene assegnato il valore nullo
         *  I Flags sono inizializzati tutti con FALSE
         *  Gli oggetti inizializzati ma vuoti
         *  il numero di threads viene inizializzato invece con il valore passato da terminale.
         */
        HDP_MCMC();

        /*! \brief Distruttore di default
         */
        ~HDP_MCMC()=default;
        //METODI DI INIZIALIZZAZIONE DELLE VARIABILI

        /*! \brief Imposta il numero di cluster iniziale che desidera l'utente
         *  \param _K - Cluster iniziali
         */
        void SetK_init(unsigned int _K);

        /*! \brief Acquisisce il dataset da file e le dimensioni del dataset
         *  Chiama dei metodi di Corpus che si occupano di creare le struttre che gestiscono il dataset
         *  \param Dataset - Nome del file che contiene i dati
         *  \param MainVariable - Nome del file che contiene le informazioni su delle dimensioni
         */
        void SetDataset(const std::string& Dataset,const std::string& MainVariable);

        /*! \brief Imposta Alpha fisso, il cui valore è deciso dall'utente
         *  \param  _Alpha - Valore fisso assegnato ad alpha dall'utente
         */
        void SetAlphaFixed(double _Alpha);

        /*! \brief Imposta una prior su Alpha
         *  \param _AA - parametro di forma passato dall'utente
         *  \param _AB - parametro di rate passato dall'utente
         */
        void SetAlphaPrior(double _AA, double _AB);

         /*! \brief Imposta Gamma fisso, il cui valore è deciso dall'utente
         *  \param  _Gamma - Valore fisso assegnato a Gamma dall'utente
         */
        void SetGammaFixed(double _Gamma);

        /*! \brief Imposta una prior su Gamma
         *  \param _GA - parametro di forma passato dall'utente
         *  \param _GB - parametro di rate passato dall'utente
         */
        void SetGammaPrior(double _GA, double _GB);

        /*! \brief Imposta il valore su Lambda
         *  \param Lambda - Parametro passato dall'utente
         */
        void SetLambdaInfo(HYP Lambda);

        /*! \brief Fissa il seed, unico per l'intera esecuzione dell'algoritmo
		 *  \param Seed - seed
		 */
        void SetSeed(const unsigned long Seed);

        /*! \brief Impone il controllo dell'algoritmo con LPML
         *  \param burnin - indica da quale iterazioni inziare a calcolarlo
         */
        void Check_Model(unsigned long burnin);

        /*! \brief Estra la dimensione del parametro latente
         *  \return dimensione parametro latente
         */
        unsigned int ViewW();

        /*! \brief Estra il numero di cluster attuale
         *  \return Numero di cluster
         */
        unsigned int ViewK();

        /*! \brief Estra il numero totale di gruppi stanziato
         *  \return Numero totale di gruppi
         */
        unsigned int ViewD();

        /*! \brief Estra la dimensione del dataset
         *  \return dimensione del dataset
         */
        unsigned int ViewN();

        /*! \brief	Algoritmo globale vedi Algorithm 3 Relazione_Parisi_Perego.pdf
		 *  \param Iterations - numero massimo di iterazioni (criterio di stop)
		 *  \param Iterations_Sub - numero massimo di iterazione per il gibbs samplere che campiona i sub-topic (criterio di stop)
		 */
        void Algorithm(unsigned int Iterations, unsigned int Iterations_Sub);


    private:
        /*! \brief Scambia le vecchie proposte inserendole nelle strutture definitive con quelle nuove
         * Quando il metodo si chiude, le nuove proposte, che non sono state accettate nei passi di
         * M-H, sono distrutte. Nelle mosse globali, le modifiche per le nuove proposte per i cluster e per le etichette vengono
         * fatte direttamente nelle strutture definitive. Nel caso queste non sono accettate è necessario reinserie nelle struttire definitive
         * la situazione precedente alle proposte.
         * E' un metodo template perche' e' usato sia per la struttura che racchiude i gruppi, Corpus, sia per la struttra che
         * gestisce i clusters, Model.
         *  \param Old - Situazione antecedente al cambiamento
         *  \param New - Nuove proposte, non accettate
         */
        template <class T> void Swap(T& Old,T& New);

        /*! \brief Stampa a video un sunto delle impostazioni dell'algoritmo scelte
         */
        void Summary();

        /*! \brief Inizializza la struttura che gestisce i Clusters con il numero di cluster iniziale scelto
         */
        void SetClusters();

        /*! \brief Aggiorna gli iperparametri del parametro lantete, in base a come sono distribuiti i dati nei vari clusters
         */
        void UpdateClusterCounts();

        /*! \brief aggiorna gli \f$ m_{jk} \f$ in ogni gruppo con l'equazione 3.3 in Relazione_Parisi_Perego.pdf
         *  Dopo di che \f$ m_{.k} = \sum_j m_{jk} \f$ in modo da aggiornare i tavoli nel clusters
         */
        void UpdateTable();

        /*! \brief aggiorna i tavoli nei sub.cluster di ogni clusters, \f$ m_{jkl},  m_{jkr} \f$ in ogni gruppo
         *  con l'equazione 3.12 in Relazione_Parisi_Perego.pdf
         *  Dopo di che \f$ m_{.kh} = \sum_j m_{jkh} \f$ in modo da aggiornare i tavoli in dei sub-clusters di ogni cluster
         */
        void UpdateSubTable();

        /*! \brief Aggiorna in ogni gruppo i pesi dei clusters che compiaiono nel gruppo
         */
        void UpdateDocWeights();

        /*! \brief Aggiorna in ogni gruppo i pesi dei sub-clusters di ogni cluster che compare nel gruppo
         */
        void UpdateDocWeights_Sub();

        /*! \brief assegna le nuove etichette del cluster e del sub-cluster ad ogni dato
         */
        void UpdateAssignment_Cluster_and_Subcluster();

        /*! \brief assegna le nuove etichette del cluster, senza assegnare il sub-cluster, ad ogni dato
         * Viene utilizzato nelle mosse di Merge/Split globale per fare le nuove prposte per le nuove etichette
         */
        void UpdateAssignment_Cluster();

        /*! \brief Aggiorna i pesi globali di ogni cluser, compreso quello vuoto;
         *  Equazione di campionamento 3.4 di Relazione_Parisi_Perego.pdf
         */
        void UpdateBeta();

        /*! \brief Aggiorna i pesi globali dei sub-cluster di ogni cluser, compresi i sub-cluster del cluster vuoto;
         *  Equazione di campionamento 3.8 di Relazione_Parisi_Perego.pdf
         */
        void UpdateAllBetaSub();

         /*! \brief Aggiorna i pesi globali dei sub-cluster del cluser identificato dal suo ID
          *  Utilizzato nel metodo UpdateAllbetaSub() e nei passi di Gibbs-sampler per le proposte dei sub-topic
          *  dopo aver accetatto le mosse di M-H.
          *  Equazione di campionamento 3.8 di Relazione_Parisi_Perego.pdf
          *  \param k - id del cluster di cui si vuole aggiornare i pesi globali dei sui sub-cluster
          */
        void UpdateBetaSub(const ClusterID k);
        //usa il singolo dato
        //double Loglikelihood(const ClusterID _k);
        //void Loglikelihood_SubCluster(const ClusterID _k, pair<double,double>& LogLike);
        /*! \brief Controlla quali cluster sono vuoti e li elimina
         */
        void EmptyCluster();

        /*! \brief Controlla se uno dei due sub-cluster del cluster identificato con il suo ID, e' vuoto.
         *  \return TRUE se uno dei sub-cluster e' vuoto.
         */
        bool IsEmptySubcluster(const ClusterID _k);

        /*! \brief Calcola la matrice 8.12 con l'equazione 8.13, vedi Relazione_Parisi_Perego.pdf
         */
        void computeLogL();

        /*! \brief Calcola la quantità 3.25 della Relazione_Parisi_Perego.pdf
         */
        long double logq();

        /*! \brief Passi di Gibbs Sampler per campionare le nuove proposte per i sub-cluster
         *  dei nuovi cluster. Vengono utilizzate le equazioni da 3.8 a 3.12
         *  \param ProposedClusters - Etichette dei nuovi cluster di cui bisogna fare la proposta per i rispettivi sub-cluster
         */
        void Gibbs_SubCluster (const vector<ClusterID>& ProposedClusters);

        /*! \brief Mosse di split locale.
         *  Si propone la divisione di un cluster alla volta, in due cluster.
         *  Questa proposta è la stessa per tutti i gruppi.
         *  Se si accetta la nuova prposta si campionano si sub-cluster dei due nuovi cluster in ogni gruppo
         */
        void LocalSplit();

        /*! \brief Mosse di merge locale.
         *  Si formano delle coppie casuali dei cluster attualmente presenti, cui proporre l'unione per formare un
         *   nuovo cluster.
         *  Se si accetta la nuova proposta si campionano i sub-cluster del nuovo cluster in ogni gruppo.
         */
        void LocalMerge();

        /*! \brief Mosse di merge globale.
         *  Si formano delle una coppia casuale dei cluster attualmente presenti, cui proporre l'unione per formare un
         *   nuovo cluster.
         *  Si campionano le nuove quantità per tutti i cluster, anche per quelli che erano già presenti.
         *  Se si accetta la nuova proposta si campionano i sub-cluster di tutti i clusters.
         */
        void GlobalMerge();

        /*! \brief Mosse di split globale.
         *  Si sceglie casualmente un cluster non vuoto di cui proporre lo split e formare due nuovi clusters
         *  Si campionano le nuove quantità per tutti i cluster, anche per quelli che erano già presenti.
         *  Se si accetta la nuova proposta si campionano i sub-cluster di tutti i clusters.
         */
        void GlobalSplit();

        /*! \brief Campionamento della nuova Alpha, da usare nell'iterazione successiva dell'algoritmo
         */
        void AlphaPrior();

        /*! \brief Campionamento della nuova Gamma, da usare nell'iterazione successiva dell'algoritmo.
         */
        void GammaPrior();

        /*! \brief Aggiornamento del K corrente a fronte di aggiunta o eliminazione di clusters
         */
        void UpdateK();

        /*! \brief Aggiornamento dello storico dei K. Fatto alla fine di ogni iterazione
         */
        void UpdateAllK();

        /*! \brief Salva su file lo storico dei K  alla fine dell'algoritmo, per fare le analisi a posteriori
         */
        void SaveAllK();

        /*! \brief Salva su file lo storico degli Alpha, alla fine dell'algoritmo, per fare le analisi a posteriori
         */
        void SaveAllAlpha();

        /*! \brief Salva su file lo storico dei Gamma, alla fine dell'algoritmo, per fare le analisi a posteriori
         */
        void SaveAllGamma();

        void SaveRunTime();

        /*! \brief Salva su file i pesi globali dei clusters dopo 10000 iterazioni ogni 10 iterazioni
         */
        void SaveLastBeta(const std::string& Filename);

		/*! \brief Salva su file i pesi dei cluster specifici dei documenti dopo 10000 iterazioni ogni 10 iterazioni
		*/
		void SaveLastPi(const std::string& Filename);

        /*! \brief Salva su file i parametri lantenti di ogni cluster dopo 10000 iterazioni ogni 10 iterazioni
         */
        void SaveLastTheta(const std::string& Filename);

        /*! \brief Salva su file le quantià necessarie per il calcolo dell' LPML, scartanto le prime burnin iterazioni
         */
        void LPML();

        /*! \brief Salva su file le etichette dei dati, dopo 10000 iterazioni ogni 10 iterazioni
         */
        void SaveLabels(const std::string& Filename);

		//void ReadPi();

	};

/*
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::ReadPi(){

	std::ifstream LastPi("./cpp_results/LastPi.bin",std::ios::binary);

    if(LastPi.fail()){
		std::cerr<<"Doesn't open file: LastPi"<<std::endl;
		exit(1);
    }

    //while(check){

    unsigned int ITER_Pi;
    LastPi.read(reinterpret_cast<char*>(&ITER_Pi), sizeof(unsigned int));

	unsigned int iter = ITER_Pi;


    for(size_t it = 0; it<MaxIt - ITER_Pi; it++){

        unsigned int K;
        LastPi.read(reinterpret_cast<char*>(&K), sizeof(unsigned int));

		std::cout<<"iter "<<iter<<std::endl;

		for(size_t d=0; d<D; ++d){
			double *tmp_pi = new double [K+1];
			LastPi.read(reinterpret_cast<char*>(tmp_pi), (K+1)*sizeof(double));

			std::vector<double> TempPi(tmp_pi, tmp_pi + K+1);
			delete [] tmp_pi;

			std::cout<<"d "<<d<<std::endl;

			for(auto i: TempPi)
				std::cout<<i<<" ";
			std::cout<<std::endl;
		}

		++iter;

    } //chiudo for

		//check=false;
    //} // chiudo while

    LastPi.close();

}
*/

//Costrutture
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
HDP_MCMC<MODEL,DOCUMENT,DIM>::HDP_MCMC(){

    //devo impostare tutte le informazioni dei Falgs a 0 in quanto non ho ancora acquisito nulla
    std::get<0>(Flags)=false;
    std::get<1>(Flags)=false;
    std::get<2>(Flags)=false;
    std::get<3>(Flags)=false;
    std::get<4>(Flags)=false;
    std::get<5>(Flags)=false;
	std::get<6>(Flags)=false;
    std::get<7>(Flags)=false;
    std::get<8>(Flags)=false; //LPML&Perplexity

	//nr di documenti
	D = 0;;
	//nr di parole
	N = 0;
	//Dimensione dell'iperparametro della distribuzione del parametro latente
	W = 0;
	//nr iterazione corrente
	It = 0;
	//nr massimo iterazioni (criterio stop mcmc)
	MaxIt = 0;
	//nr massimo di terazioni per criterio campionamento SubTopic, duarnte le mosse di split
	MaxIt_SubCluster = 0;
	//iperparametro massa DP globale (nel caso gamma fisso)
	Gamma = 0.0;
	AllGamma.clear();
	//iperparametri nel caso ho priori su Gamma
	GA = 0.0;
	GB = 0.0;
	//iperparametri per campionare nei documenti
	Alpha = 0.0;
	AllAlpha.clear();
	//iperparamtri nel caso ho prior su Alpha
	AA = 0.0;
	AB = 0.0;
	//storico dei K
	AllK.clear();

	LogStirlingNumbers.clear();

	m = 0;

	CPO.clear();
	burnin_CPO = 0;
    // Valore di default numero massimo di thread (seriale)
	OMP_NUM_THREADS=1;
	// Lettura numero massimo di thread
	char* EnvironmentVariable;
	EnvironmentVariable = getenv("OMP_NUM_THREADS");
	if (EnvironmentVariable!=NULL){
		std::stringstream ss;
		ss << EnvironmentVariable;
		ss >> OMP_NUM_THREADS;
	}


}

//Fissa i K iniziali
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetK_init(unsigned int _K){

	K = _K;
	std::get<6>(Flags) = true;  // l'utente vorrebbe impostare K


}

//Aggiorna i K correnti
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateK(){

    K=Model.ViewK();
}

//Aggiorna lo storico dei K
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateAllK(){

    this->UpdateK();
    AllK.push_back(K);
}

//Fissa il seme
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetSeed(const unsigned long Seed){

    if(It != 0){
        std::cerr<<"Warning: Seed can be set only at the beggining of the algorithm"<<std::endl;
        return;
    }

    if(get<7>(Flags))
        std::cerr<<"Warning: Overwritng of the seed in the algorithm"<<std::endl;

    Gen.fixedSeed(Seed);
	get<7>(Flags)=true;


	return;


}

//Mostra la dimensione dell'iperparametro della distribuzione del parametro latente
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
unsigned int HDP_MCMC<MODEL,DOCUMENT,DIM>::ViewW(){

    return W;
}

//Mostra il K corrente
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
unsigned int HDP_MCMC<MODEL,DOCUMENT,DIM>::ViewK(){

    return K;
}

//Mostra quanti gruppi abbiamo
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
unsigned int HDP_MCMC<MODEL,DOCUMENT,DIM>::ViewD(){

    return D;
}

//Mostra la dimensione del dataset
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
unsigned int HDP_MCMC<MODEL,DOCUMENT,DIM>::ViewN(){

    return N;

}

//Acquisizione Dataset
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetDataset(const std::string& Dataset,const std::string& MainVariable){

    if(It!=0){
        std::cerr<<"ERROR Dataset can be set only before the running of Algorithm"<<std::endl;
        exit(1);
    }

    if(std::get<0>(Flags)){
        std::cerr<<"WARNING overwriting Dataset"<<std::endl;
    }

    std::ifstream file1 (MainVariable);
    vector<unsigned int> Variable;

    if(file1.fail()){
        std::cerr<<"ERROR in Variables acquisition"<<std::endl;
        std::cerr<<"I can't open: "<<MainVariable<<std::endl;
        exit(1);
    }

    std::string ss1;
    unsigned int temp;

    while(std::getline(file1,ss1)){
        std::istringstream SSTR(ss1);
        SSTR>>temp;
        Variable.push_back(temp);
    }

    if(Variable.size()!=3){
        std::cerr<<"ERROR in acquisition Variables"<<std::endl;
        std::cerr<<"Wrong number of variables"<<std::endl;
        exit(1);
    }

    D=Variable[0];
    W=Variable[1];
    N=Variable[2];

    this-> SetClusters();

	std::ifstream file (Dataset);

    if(file.fail()){
        std::cerr<<"Error in Dataset acquisition"<<std::endl;
        std::cerr<<"I cannot open "<<Dataset<<std::endl;
        exit(1);
    }


	vector<unsigned int> AllNj;   // per controllare il nr totale di parole inserite
	AllNj.reserve(D);
	corpus.reserve(D);

	DOCUMENT<DIM> Doc(Alpha);

	for(typename Corpus::size_type d=0; d<D; ++d)
		corpus.push_back(Doc);

	std::string ss;
	unsigned int DocId;

	while(std::getline(file,ss)){
        std::istringstream SSTR(ss);
        SSTR>>DocId;
		corpus[DocId-1].SetDataset(SSTR);
	}

	unsigned int Nj;

	for(typename Corpus::size_type d=0; d<D; ++d){
		Nj = corpus[d].SortData(K,Gen);
		AllNj.push_back(Nj);
	}

	//controllo sul nr di parole inserite
	unsigned int sum=0;

	for (auto i: AllNj)
	     sum+=i;

	if(sum!=N){
	  std::cerr<<" Total number of data inserted doesn't match N"<<std::endl;
	  exit(1);
	}

    unsigned int Nmax=0;
	auto it_max=std::max_element(AllNj.cbegin(),AllNj.cend());
	Nmax= *(it_max);
	unsigned int dim= (Nmax+1)*(Nmax+2)/2;
	LogStirlingNumbers.resize(dim);

	ComputeLogStirlingNumbers (Nmax, LogStirlingNumbers);

	std::get<0>(Flags)=true;   //ho impostato dataset

	this-> UpdateClusterCounts();

	this-> EmptyCluster();


}

//Fissa il numero di cluster iniziali anche nella struttura Model
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetClusters(){

    if(std::get<6>(Flags)){
        Beta_empty = 1.0/(K+1);
        Model.SetInitialClusters(K);
    }
    else {
        Gen.setNumThreads(1);
        K = Gen.runifdiscrete(10) + 1;
        Beta_empty = 1.0/(K+1);
        Model.SetInitialClusters(K);
    }


}

//Stabilisce se bisogna utilizzare LPML
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::Check_Model(unsigned long burnin){

    std::get<8>(Flags) = true;
    burnin_CPO = burnin;

}

//Metodo che permette do effettuare uno scambio
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
 template<class T> void HDP_MCMC<MODEL,DOCUMENT,DIM>::Swap(T& Old, T& New){

    T Temp(std::move(Old));
    Old = std::move(New);
    New = std::move(Temp);

}

//Stampa a video le informazioni sull'esecuzione dell'algoritmo
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::Summary(){

    std::cout<<"**********************************************"<<std::endl;
    std::cout<<"          Informazioni riassuntive:           "<<std::endl;
    std::cout<<"**********************************************"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Numero di gruppi: "<<D<<std::endl;
    std::cout<<"Dimensione parametro latente: "<<W<<std::endl;
    std::cout<<"Dimensione Dataset: "<<N<<std::endl;
    std::cout<<std::endl;

    if(std::get<1>(Flags) && std::get<2>(Flags))
        std::cout<<"Alpha ~ Gamma("<<AA<<" , "<<AB<<")"<<std::endl;

    if((std::get<1>(Flags) && !std::get<2>(Flags)) || (!std::get<1>(Flags)))
        std::cout<<"Alpha = "<<Alpha<<std::endl;

    if(std::get<3>(Flags) && std::get<4>(Flags))
        std::cout<<"Gamma ~ Gamma("<<GA<<" , "<<GB<<")"<<std::endl;

    if((std::get<3>(Flags) && !std::get<4>(Flags)) || (!std::get<3>(Flags)))
        std::cout<<"Gamma = "<<Gamma<<std::endl;

    Model.PrintLambdaInfo();
    std::cout<<std::endl;

    std::cout<<"Iterazioni: "<<MaxIt<<std::endl;
    std::cout<<"Iterazioni sub-Cluster: "<<MaxIt_SubCluster<<std::endl;
    std::cout<<"Clusters iniziali: "<<K<<std::endl;
    std::cout<<"OMP_NUM_THREADS: "<<OMP_NUM_THREADS<<std::endl;
	std::cout<<std::endl;

}

//Acquisisce l'informazione riguardante l'iperparametro Lamda
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetLambdaInfo(HYP Lambda){

    if(It!=0){
        std::cerr<<"Information on Lambda can be set only at the begging of the Algoritm"<<std::endl;
        return;
    }

    if(std::get<5>(Flags))
        std::cerr<<"WARNING Overwriting information on Lambda. It has been already inserted"<<std::endl;

    Model.SetHyperparameter(Lambda);

    std::get<5>(Flags)=true;


}

//Acquisisce le informazioni di Alpha nel caso in cui si decida di avrlo fisso
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetAlphaFixed(double _Alpha){

    if(It!=0){
        std::cerr<<"Information on Alpha can be set only at the begging of the Algoritm"<<std::endl;
        return;
    }

    if(std::get<1>(Flags))
        std::cerr<<"WARNING Overwriting information on Alpha. It has been already inserted"<<std::endl;

    if(std::get<1>(Flags) && std::get<2>(Flags))
        std::cerr<<"WARNING overwriting Alpha prior, with fixed value"<<std::endl;


    if(_Alpha<0){
        std::cerr<<"ERROR Alpha has to be Positive"<<std::endl;
        exit(1);
    }

    Alpha= _Alpha;
    AA=0;
    AB=0;

    std::get<1>(Flags)=true;
    std::get<2>(Flags)=false;


}

//Acquisisce le informazioni della Prior su Alpha
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetAlphaPrior(double _AA, double _AB){

    if(It!=0){
        std::cerr<<"Information on Alpha can be set only at the begging of the Algoritm"<<std::endl;
        return;
    }

    if(std::get<1>(Flags))
        std::cerr<<"WARNING Overwriting information on Alpha. It has been already inserted"<<std::endl;

    if(std::get<1>(Flags) && !std::get<2>(Flags))
        std::cerr<<"WARNING overwriting fixed Alpha , with Prior"<<std::endl;

    if(_AA<=0){
        std::cerr<<"ERROR Generated in SetHyperameter"<<std::endl;
        std::cerr<<"Alpha shape hyperparameter must be a positive number." <<std::endl;
        exit(1);
    }

    if(_AB<=0){
        std::cerr<<"ERROR Generated in SetHyperameter"<<std::endl;
        std::cerr<<"Alpha rate hyperparameter must be a positive number." <<std::endl;
        exit(1);
    }

    AA=_AA;
    AB=_AB;
    Alpha=_AA/_AB;

    std::get<1>(Flags)=true;
    std::get<2>(Flags)=true;


}

//Aquisice informazioni sulla prior di Gamma nel caso in cui deidiamo che debba essere fissa
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetGammaFixed(double _Gamma){

    if(It!=0){
        std::cerr<<"Information on Gamma can be set only at the begging of the Algoritm"<<std::endl;
        return;
    }

    if(std::get<3>(Flags))
        std::cerr<<"WARNING Overwriting information on Gamma. It has been already inserted"<<std::endl;

    if(std::get<3>(Flags) && std::get<4>(Flags))
        std::cerr<<"WARNING overwriting Gamma prior, with fixed value"<<std::endl;

    if(_Gamma<0){
        std::cerr<<"ERROR Gamma has to be Positive"<<std::endl;
	exit(1);
    }

    Gamma= _Gamma;
    GA=0;
    GB=0;

    std::get<3>(Flags)=true;
    std::get<4>(Flags)=false;


}

//Acquisisce le informazioni nel caso in cui la Gamma non sia fissa
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SetGammaPrior(double _GA, double _GB){

    if(It!=0){
        std::cerr<<"Information on Gamma can be set only at the begging of the Algoritm"<<std::endl;
        return;
    }

    if(std::get<3>(Flags))
        std::cerr<<"WARNING Overwriting information on Gamma. It has been already inserted"<<std::endl;


    if(std::get<3>(Flags) && !std::get<4>(Flags))
        std::cerr<<"WARNING overwriting fixed Gamma , with prior"<<std::endl;

    if(_GA<=0){
        std::cerr<<"ERROR Generated in SetHyperameter"<<std::endl;
        std::cerr<<"Gamma shape hyperparameter must be a positive number." <<std::endl;
        exit(1);
    }

    if(_GB<=0){
        std::cerr<<"ERROR Generated in SetHyperameter"<<std::endl;
        std::cerr<<"Gamma rate hyperparameter must be a positive number." <<std::endl;
        exit(1);
    }

    GA=_GA;
    GB=_GB;
    Gamma=_GA/_GB ;

    std::get<3>(Flags)=true;
    std::get<4>(Flags)=true;


}

//Campiona le Gamma da utilizzare nell'iterazione successiva
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::GammaPrior(){

	Gen.setNumThreads(1);

	if(!(m>0)){
	  std::cerr<< "error in Gamma prior , parametri Beta >0"<<std::endl;
	  exit(1);
	}

	double eta = Gen.rbeta(Gamma+1.0,m);           //variabile ausiliaria
	double R = (Gamma + K -1.0)/(m*(GB-log(eta)));
    double pi_eta = R / (1.0+R);

	double u = Gen.runif();

	if( !((GA+K)>0.0 ) || !(1.0/(GB-log(eta))>0) || !((GA+K-1.0)>0 )){
	  std::cerr<< "error in Gamma prior , parametri Gamma >0"<<std::endl;
	  exit(1);
	}

	if(u<pi_eta){
	  Gamma = Gen.rgamma(GA+K,1.0/(GB-log(eta)));
	}
	else{
	  Gamma = Gen.rgamma(GA+K-1.0,1.0/(GB-log(eta)));
	}

	AllGamma.push_back(Gamma);


}

//Campiona le Alpha da utilizzare nell'iterazione successiva
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>    // Teh code, randutils.c, randconparam
void HDP_MCMC<MODEL,DOCUMENT,DIM>::AlphaPrior(){

	vector<unsigned int> DatasInGroups(D);
    vector<double> wj(D);
    vector<unsigned int> sj(D);

	Gen.setNumThreads(OMP_NUM_THREADS);

	typename Corpus::size_type d;
	unsigned int sum_sj = 0;
	double sum_log_wj = 0;

	for(size_t it = 0; it<20; ++it){

	#pragma omp parallel default(none) private(d) shared(DatasInGroups,sum_sj,sum_log_wj,sj,wj,cerr)
	{
	#pragma omp for schedule(static,1) reduction( +: sum_sj,sum_log_wj)
    for(d=0; d<D; d++){
        DatasInGroups[d] = corpus[d].ViewNj();

		if( !( DatasInGroups[d] > 0)){
		  std::cerr<< "error in Alpha Prior, beta parameter > 0"<<std::endl;
		  exit(1);
		}

		wj[d] = Gen.rbeta(Alpha+1.0,DatasInGroups[d]) ;

		if( !((DatasInGroups[d]/(Alpha+ DatasInGroups[d])) > 0.0) || !((DatasInGroups[d]/(Alpha+DatasInGroups[d])) <= 1.0) ){
		  std::cerr<< "error in Alpha Prior, bernoulli parameter"<<std::endl;
		  exit(1);
		}

        sj[d] = Gen.rbernoulli( DatasInGroups[d]/(Alpha + DatasInGroups[d] ));
		sum_log_wj += std::log(wj[d]);
		sum_sj += sj[d];

    }
	}

	if( !( (AA + m - sum_sj) > 0) || !( (1.0/( AB - sum_log_wj)) > 0)){
	  std::cerr<< " error in alpha prior, parametri gamma "<<std::endl;
	  exit(1);
	}

    Alpha = Gen.rgamma( AA + m - sum_sj, 1.0/( AB - sum_log_wj));

	sum_log_wj = 0.0;
	sum_sj = 0;
	}

    AllAlpha.push_back(Alpha);


}

//Algoritmo
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::Algorithm(unsigned int Iterations, unsigned int Iterations_Sub){

	if(Iterations >= It){
        MaxIt = Iterations;
        AllK.reserve(Iterations + 1);

        if(std::get<1>(Flags) && std::get<2>(Flags))
            AllAlpha.reserve(Iterations + 1);

        if(std::get<3>(Flags) && std::get<4>(Flags))
            AllGamma.reserve(Iterations + 1);
    }
    else{
        std::cerr <<"ERROR GENERATED in Algorithm" <<std::endl;
		std::cerr <<"Algorithm isn't running (low number of iterations)." <<std::endl;
		exit(1);
    }

    if(!std::get<0>(Flags)){
        std::cerr<<"ERROR Algorithm isn't running"<<std::endl;
        std::cerr<<"Dataset wasn't inserted"<<std::endl;
        exit(1);
    }

    if(!std::get<1>(Flags))
        Alpha = 1.0;

    if(!std::get<3>(Flags))
        Gamma = 1.0;

   if(!std::get<5>(Flags))
        Model.DefaultHyperparameter(W);

    MaxIt_SubCluster = Iterations_Sub;

	AllK.push_back(K);

	if(std::get<1>(Flags) && std::get<2>(Flags)){
	   AllAlpha.push_back(Alpha);
	}

	if(std::get<3>(Flags) && std::get<4>(Flags)){
	   AllGamma.push_back(Gamma);
	}

	if(std::get<8>(Flags)){
        CPO.assign(N,0.0);
	}

    this->Summary();

    //////////////////////////////////////////
    //////////////FILE NAME///////////////////
    //////////////////////////////////////////
    std::string PathTheta = "./cpp_results/Theta";
    std::string PathPi = "./cpp_results/Pi";
    std::string PathLabels = "./cpp_results/Labels";
    std::string PathBeta = "./cpp_results/Beta";
    std::string extention = ".bin";
    std::string ThetaFile = "./cpp_results/Theta.bin";
    std::string PiFile = "./cpp_results/Pi.bin";
    std::string LabelsFile = "./cpp_results/Labels.bin";
    std::string BetaFile = "./cpp_results/Beta.bin";
    ///////////////////////////////////////////

    for(It =0; It < MaxIt; It++ ){

//        if (It%100 == 0)
            std::cout<<"I'm doing iteration "<<It<<"..."<<std::endl;
std::cout<<K<<std::endl;
        this-> UpdateDocWeights();

        this-> UpdateDocWeights_Sub();

        Model.UpdateThetaCluster(Gen);

        Model.UpdateThetaSubCluster(Gen);

        this-> UpdateAssignment_Cluster_and_Subcluster();

        this-> EmptyCluster();

        this -> computeLogL();

        this-> LocalMerge();

        this-> EmptyCluster();

        this-> LocalSplit();

        this-> EmptyCluster();

        this-> GlobalMerge();

        this-> EmptyCluster();


        this-> GlobalSplit();

        this-> EmptyCluster();

        this-> UpdateTable();

        this-> UpdateSubTable();

        this-> UpdateBeta();

        this-> UpdateAllBetaSub();


        //if(MaxIt - It <= 100){
        if(It >=1000){
        //if(It%10==0){
           // std::cout<<"Sto salvando l'iterazione "<<It<<std::endl;
            if(It%500==0){
                std::cout<<"Sto cambiando nome ai file"<<std::endl;
                std::stringstream add;

                add<<PathBeta<<It<<extention;
                BetaFile.clear();
                BetaFile = add.str();
                std::cout<<"Sto salvando su file "<<BetaFile<<std::endl;
                this->SaveLastBeta(BetaFile);

                add.str(std::string());

                add<<PathPi<<It<<extention;
                PiFile.clear();
                PiFile = add.str();
                std::cout<<"Sto salvando su file "<<PiFile<<std::endl;

                this->SaveLastPi(PiFile);

                add.str(std::string());

                add<<PathTheta<<It<<extention;
                ThetaFile.clear();
                ThetaFile = add.str();
			                std::cout<<"Sto salvando su file "<<ThetaFile<<std::endl;

                this->SaveLastTheta(ThetaFile);

                add.str(std::string());

                add<<PathLabels<<It<<extention;
                LabelsFile.clear();
                LabelsFile = add.str();
                                std::cout<<"Sto salvando su file "<<LabelsFile<<std::endl;

                this->SaveLabels(LabelsFile);

                add.str(std::string());
			}

			else{
                                std::cout<<"Sto salvando su file "<<BetaFile<<std::endl;
                                std::cout<<"Sto salvando su file "<<PiFile<<std::endl;
                                std::cout<<"Sto salvando su file "<<LabelsFile<<std::endl;
                                std::cout<<"Sto salvando su file "<<ThetaFile<<std::endl;

                this -> SaveLastBeta(BetaFile);
                this ->SaveLastPi(PiFile);
                this ->SaveLabels(LabelsFile);
                this ->SaveLastTheta(ThetaFile);
			}

        }

        //if(It > 5000)
        //    this->SaveLastTheta();

		this -> UpdateAllK();

		if(std::get<1>(Flags) && std::get<2>(Flags))
            this->AlphaPrior();

        if(std::get<3>(Flags) && std::get<4>(Flags))
            this->GammaPrior();


		//if(N<10001 && MaxIt - It <= 10000)
		//if(It > 5000)
			//this -> SaveLabels();

        if(std::get<8>(Flags) && It >= burnin_CPO)
            this -> LPML();
        std::cout<<"Sto salvando RunTime"<<std::endl;
        this -> SaveRunTime();
    }

	//this->SaveAllK();

	//if(std::get<1>(Flags) && std::get<2>(Flags))
	  // this->SaveAllAlpha();

	//if(std::get<3>(Flags) && std::get<4>(Flags))
	  // this->SaveAllGamma();

	//this-> ReadPi();

}

//Aggiorna i conteggi di ogni cluster,quanti elementi sono caduti in ogni cluster
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateClusterCounts(){

    omp_set_num_threads(OMP_NUM_THREADS);

	#pragma omp parallel
	{

        STAT counts4cleft(W);
        STAT counts4cright(W);

        #pragma omp for schedule(static,1)
        for (size_t k=0; k<K; ++k){

            Model[k].ResetStatistics(W);
            Model[k].ResetStatisticsLeft(W);
            Model[k].ResetStatisticsRight(W);

            for(typename Corpus::size_type j=0; j<D; ++j){

                corpus[j].ViewCounts4c(k,counts4cleft,counts4cright);

                Model[k].UpdateStatistics(counts4cleft,counts4cright);
            }

        }

	}

}

//Campiona i tavoli di tutti i cluster
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateTable(){

    unsigned int sum_mk;
    vector<double> _Beta;

    Model.ViewBeta(_Beta);
    Gen.setNumThreads(OMP_NUM_THREADS);

    #pragma omp parallel for schedule(static,1)
	for(typename Corpus::size_type d=0; d<D; d++){
		corpus[d].UpdateLocalTable( LogStirlingNumbers , _Beta, Gen );
	}

    m = 0;

    for(size_t i = 0; i < K; i++){
        sum_mk = 0;

        for(typename Corpus::size_type d=0;d<D;d++)
            sum_mk += corpus[d].ViewNumTableID(i);

        Model[i].SetGlobalTable(sum_mk);

		m += sum_mk; // nr totale di tavoli
    }

}

//Campiona i tavoli dei sub-cluster di ogni cluster
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateSubTable(){

    pair<unsigned int, unsigned int> sum_mk;
    vector<double> _BetaL;
    vector<double> _BetaR;

    Model.ViewBetaLeft(_BetaL);
    Model.ViewBetaRight(_BetaR);

    Gen.setNumThreads(OMP_NUM_THREADS);

    #pragma omp parallel for schedule(static,1)
    for(typename Corpus::size_type d=0; d<D; d++){
        corpus[d].UpdateAllLocalTableSub(LogStirlingNumbers, _BetaL, _BetaR, Gen);
    }

	for(size_t i = 0; i < K; i++){

        sum_mk=std::make_pair(0,0);

        for(auto it = corpus.begin(); it != corpus.end(); it++){
            sum_mk.first += (*it).ViewNumTableLeftID(i);
            sum_mk.second += (*it).ViewNumTableRightID(i);
        }

        Model[i].SetGlobalTableLeft(sum_mk.first);
        Model[i].SetGlobalTableRight(sum_mk.second);
    }

}

//Campiona i pesi dei cluster in ogni documento
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateDocWeights(){

    vector<double> AllBeta;
    AllBeta.reserve(K+1);

    Model.ViewBeta(AllBeta);

    AllBeta.push_back(Beta_empty);

    Gen.setNumThreads(OMP_NUM_THREADS);

    #pragma omp parallel for default(none) shared(AllBeta) schedule(static,1)
    for

	(typename Corpus::size_type d=0; d<D;++d){

		corpus[d].UpdatePi(AllBeta,Gen);

     }

}

//Campiona i pesi de sub-cluster in ogni cluster
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateDocWeights_Sub(){

    vector<double> BetaLeft;
    BetaLeft.reserve(K);
    vector<double> BetaRight;
    BetaRight.reserve(K);

    Model.ViewBetaLeft(BetaLeft);
    Model.ViewBetaRight(BetaRight);

    Gen.setNumThreads(OMP_NUM_THREADS);

    #pragma omp parallel for default(none) shared(BetaLeft,BetaRight) schedule(static,1)
    for(typename Corpus::size_type d=0 ; d< D; ++d){
        corpus[d].UpdateAllPiSub(BetaLeft,BetaRight,Gen);
    }

}

//Assegna le nuove etichette del cluster, senza assegnare il sub-cluster, ad ogni dato
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateAssignment_Cluster(){

    Gen.setNumThreads(OMP_NUM_THREADS);

    #pragma omp parallel default(none)
    {
        vector<POINT> VettId;
        vector<double> Likelihood;
        Likelihood.resize(K);

        #pragma omp for schedule(static,1) private(VettId)
        for(typename Corpus::size_type d=0;d<D;++d){
            corpus[d].ViewData(VettId);

            for(auto VettId_it = VettId.cbegin(); VettId_it != VettId.cend(); ++ VettId_it){

                for(vector<double>::size_type k = 0; k<K; k++)
                    Likelihood[k] = std::exp(Model.Loglikelihood((*(VettId_it)), k));

                corpus[d].UpdateZeta(Likelihood,*(VettId_it),Gen);
            }

            VettId.clear();
        }
    }

    this-> UpdateClusterCounts();

}

//Assegna le nuove etichette del cluster, assegnando il sub-cluster, ad ogni dato
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateAssignment_Cluster_and_Subcluster(){

	vector<POINT> VettId;

	Gen.setNumThreads(OMP_NUM_THREADS);

    #pragma omp parallel
    {
        vector<double> Likelihood;
        vector<double> LikelihoodLeft;
        vector<double> LikelihoodRight;
        Likelihood.resize(K);
        LikelihoodLeft.resize(K);
        LikelihoodRight.resize(K);

        #pragma omp for private(VettId) schedule(static,1)
        for(typename Corpus::size_type d=0; d<D; d++){

            corpus[d].ViewData(VettId);
            // per ogni parola nel doc d devo recuperare \theta_{id,k}, \theta_{id,kl}, \theta_{id,kr} per ogni k
            for(auto VettId_it = VettId.cbegin(); VettId_it != VettId.cend(); ++ VettId_it ) {

                for(vector<double>::size_type k = 0; k<K; k++){
                    Likelihood[k] = std::exp(Model.Loglikelihood((*(VettId_it)), k));
                    LikelihoodLeft[k] = std::exp(Model.LoglikelihoodLeft((*(VettId_it)), k));
                    LikelihoodRight[k] = std::exp(Model.LoglikelihoodRight((*(VettId_it)), k));
                }

                corpus[d].UpdateZeta_and_Sub(Likelihood,LikelihoodLeft,LikelihoodRight,*(VettId_it),Gen);
            }

            VettId.clear();
        }
    }

    this-> UpdateClusterCounts();

}

//Aggiorna i pesi globali dei cluster
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateBeta(){

    Gen.setNumThreads(1);
    vector<double> params;
    params.reserve(K+1);
    vector<double> dir_sampled(K+1);		// è il vettore dei Beta1....Betak,Betak+1
    unsigned int mk=0;

    for(size_t k=0; k<K;++k){
        mk=Model[k].ViewGlobalTable();
        params.push_back(mk);
    }

    params.push_back(Gamma);

    Gen.rdirichlet(params,dir_sampled);

	for(size_t k=0; k<K; ++k){
        Model[k].SetBeta(dir_sampled[k]);
    }

    Beta_empty = dir_sampled[K];


}

//Aggiorna i pesi dei globali dei sub-cluster di tutti i cluster
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateAllBetaSub(){

    Gen.setNumThreads(OMP_NUM_THREADS);

	#pragma omp parallel for schedule(static,1)
    for(size_t k=0; k<K; ++k)
        this->UpdateBetaSub(k);

}

//Aggiorna i pesi globale dei sub-cluster di un cluster, identificato con il suo ID
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::UpdateBetaSub(const ClusterID k){

    unsigned int mkl=0;
    unsigned int mkr=0;
    vector<double> params;
    vector<double> dir_sampled;
    params.reserve(2);
    dir_sampled.resize(2);

    mkl= Model[k].ViewGlobalTableLeft();

    params.push_back(mkl + Gamma);

    mkr= Model[k].ViewGlobalTableRight();

    params.push_back(mkr + Gamma);

    Gen.rdirichlet(params,dir_sampled);

    Model[k].SetBetaLeft(dir_sampled[0]);
    Model[k].SetBetaRight(dir_sampled[1]);


}

//Conrolla se i sub-cluster del cluster identificato dal suo id, sono vuoti
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
bool HDP_MCMC<MODEL,DOCUMENT,DIM>::IsEmptySubcluster(const unsigned int _k){

    unsigned int sum_dx = 0;
    unsigned int sum_sx = 0;

    for (auto it= corpus.begin(); it != corpus.end(); ++it){
        sum_sx += it-> CheckLeftSubcluster(_k);
        sum_dx += it-> CheckRightSubcluster(_k);
    }

    if ((sum_dx== 0) || (sum_sx==0))
        return true;
    else
        return false;


}

//Controlla se  ci sono cluater vuoti e li elimina
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::EmptyCluster(){

    vector<unsigned int> KToRemove;

    for(size_t k=0; k<K; k++){

        if(Model[k].IsEmpty())
            KToRemove.push_back(k);
    }

    if(KToRemove.empty())
        return;

	Model.RemoveClusters(KToRemove);

	for(typename Corpus::size_type d=0; d<D; d++)
	   corpus[d].RemoveCluster(KToRemove);

    this->UpdateK();


}

//Gibbs per campionare i sub-cluster dei nuovi cluster
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::Gibbs_SubCluster ( const vector<ClusterID>& ProposedClusters){

    BETA _Beta;
    vector<double> LikelihoodLeft(1);
    vector<double> LikelihoodRight(1);
    unsigned int K_proposed = ProposedClusters.size();
    NUMTABLE _m;

    vector<ClusterID>::size_type k;

    vector<pair<POINT,unsigned int>> nidjk;
    vector<ClusterID>::size_type ii;
    typename vector<pair<POINT,unsigned int>>::size_type j;

    unsigned int maleft;
    unsigned int maright;

    for(unsigned long i=0; i < MaxIt_SubCluster; i++){

        Gen.setNumThreads(OMP_NUM_THREADS);

        #pragma omp parallel for private(k) schedule(static,1)
        for(k=0; k<K_proposed; ++k)
            Model.UpdateOneThetaSubCluster(ProposedClusters[k],Gen);

        _Beta.Left.reserve(K_proposed);
        _Beta.Right.reserve(K_proposed);

        for(auto i: ProposedClusters){
            _Beta.Left.push_back(Model[i].ViewBetaLeft());
            _Beta.Right.push_back(Model[i].ViewBetaRight());
        }

        #pragma omp parallel for default(none) private(k) shared(_Beta,ProposedClusters,K_proposed) schedule(static,1)
        for(typename Corpus::size_type d = 0; d < D; d++ ){

            for(k=0; k<K_proposed; ++k)
                corpus[d].UpdatePiSub( _Beta.Left[k],_Beta.Right[k],ProposedClusters[k],Gen);  //ora updatePiSub aggiorna pil e pir di un cluster alla volta

        }

        nidjk.clear();

        #pragma omp parallel default(none) private(ii,j) firstprivate(LikelihoodLeft,LikelihoodRight,nidjk) shared(K_proposed,ProposedClusters)
        {
            #pragma omp for schedule(static,1)
            for(typename Corpus::size_type d=0; d<D;++d){
                for(ii=0; ii< K_proposed;++ii){

                    corpus[d].ViewIdCounts(nidjk,ProposedClusters[ii]);
                    corpus[d].ResetDataCountSub(ProposedClusters[ii]);

                    for(j=0; j<nidjk.size();++j){
                        LikelihoodLeft[0]= std::exp(Model.LoglikelihoodLeft( nidjk[j].first, ProposedClusters[ii]));
                        LikelihoodRight[0]=std::exp(Model.LoglikelihoodRight( nidjk[j].first, ProposedClusters[ii]));

                        corpus[d].UpdateZetaSub(LikelihoodLeft,LikelihoodRight,nidjk[j].first,nidjk[j].second, ProposedClusters[ii],Gen);
                    }

                    nidjk.clear();
                }
            }
        }


        #pragma omp parallel
        {

            STAT counts4cleft(W);
            STAT counts4cright(W);

            #pragma omp for schedule(static,1)
            for (size_t k=0; k<K_proposed; ++k){

                Model[ProposedClusters[k]].ResetStatistics(W);
                Model[ProposedClusters[k]].ResetStatisticsLeft(W);
                Model[ProposedClusters[k]].ResetStatisticsRight(W);

                for(typename Corpus::size_type j=0; j<D; ++j){

                    corpus[j].ViewCounts4c(ProposedClusters[k],counts4cleft,counts4cright);

                    Model[ProposedClusters[k]].UpdateStatistics(counts4cleft,counts4cright);
                }

            }

        }

        #pragma omp parallel for default(none) shared(_Beta,ProposedClusters,K_proposed) private(k) schedule(static,1)
        for(typename Corpus::size_type d=0 ;d<D ; ++d){
            for(k=0; k< K_proposed; ++k){
                corpus[d].UpdateLocalTableSub_OneCluster(LogStirlingNumbers, _Beta.Left[k], _Beta.Right[k], ProposedClusters[k],Gen);
            }
        }


        #pragma omp parallel for default(none) private(maleft,maright),shared(ProposedClusters,K_proposed) schedule(static,1)
        for(k=0;k<K_proposed; ++k){
            maleft = 0;
            maright = 0;

            for(typename Corpus::size_type d=0; d<D;++d){
                maleft += corpus[d].ViewNumTableLeftID(ProposedClusters[k]);
                maright += corpus[d].ViewNumTableRightID(ProposedClusters[k]);
            }

            Model[ProposedClusters[k]].SetGlobalTableLeft(maleft);
            Model[ProposedClusters[k]].SetGlobalTableRight(maright);
        }

        #pragma omp parallel for default(none) private(k) shared(ProposedClusters,K_proposed) schedule(static,1)
        for(k=0;k<K_proposed;++k)
            this -> UpdateBetaSub(ProposedClusters[k]);


        _Beta.Left.clear();
        _Beta.Right.clear();

    }

}

//Merge Locale
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT, unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::LocalMerge(){

    if(K==1)
        return;

	vector<ClusterID> CurrentClustersUnsorted;
	vector<ClusterID> CurrentClusters;
	CurrentClusters.reserve(K);
	unsigned int guess;

    Model.ViewKey(CurrentClustersUnsorted);

	Gen.setNumThreads(1);

	for(size_t i=K; i>0; --i){
	   guess = Gen.runifdiscrete(i);
	   CurrentClusters.push_back(CurrentClustersUnsorted[guess]);
	   CurrentClustersUnsorted.erase(CurrentClustersUnsorted.begin()+guess);
	}

    unsigned int NrLocalMerge =0;

    if (K%2==0)
        NrLocalMerge= K/2;
    else
        NrLocalMerge= (K-1)/2;

    double logQ_split = 0.0;
    double logQ_merge = std::log(2) + std::log(NrLocalMerge)- std::log(K)- std::log(K-1);
    vector<ClusterID> ClustersToBeRemoved;
    vector<ClusterID> ProposedClusters;
    double u = 0.0;
    ClusterID b=0;
    ClusterID c=0;
    ClusterID NewClusterID = K;

    BETA Beta;
    Beta.b_c.reserve(2);
    Beta.b_c.push_back(0.0);
    Beta.b_c.push_back(0.0);

    NJK nj;

    PI pi;
    pi.b_c.reserve(2);
    pi.b_c.push_back(0.0);
    pi.b_c.push_back(0.0);

    NUMTABLE m;
    m.Tilde_b_c.reserve(2);
    m.Tilde_b_c.push_back(0);
    m.Tilde_b_c.push_back(0);

    DATACOUNT wc;
    C counts;

    CLUSTER Cluster;
    double LogLikelihood_cap= 0.0;
    double LogLikelihood_bc= 0.0;
    long double LogQ = 0.0;
    long double LogH = 0.0;


    for(size_t i=0; i< (NrLocalMerge*2)-1; ++(++i)){

        b=CurrentClusters[i];
        c=CurrentClusters[i+1];

        Model.AddOneCluster(NewClusterID);

        Beta.b_c[0] = Model[b].ViewBeta();
        Beta.b_c[1] = Model[c].ViewBeta();
        Beta.a = Beta.b_c[0] + Beta.b_c[1];

        Model[NewClusterID].SetBeta(Beta.a);

		Model[b].ViewStatisticsLeft(counts.b_left);
		Model[b].ViewStatisticsRight(counts.b_right);

		Model[NewClusterID].UpdateStatistics(counts.b_left,counts.b_right);

		Model[c].ViewStatisticsLeft(counts.c_left);
		Model[c].ViewStatisticsRight(counts.c_right);

		Model[NewClusterID].UpdateStatistics(counts.c_left,counts.c_right);

		Gen.setNumThreads(1);

        Model.UpdateOneThetaCluster(NewClusterID,Gen);

        typename Corpus::size_type d;

        omp_set_num_threads(OMP_NUM_THREADS);

        #pragma omp parallel default(none) shared(b,c) firstprivate(pi) private(d,Cluster,nj)
        {
            #pragma omp for schedule(static,1)
            for(d =0; d<D; ++d){

                pi.b_c[0] = corpus[d].ViewPiID(b);
                pi.b_c[1] = corpus[d].ViewPiID(c);
                pi.a = pi.b_c[0] + pi.b_c[1];

                corpus[d].ViewCluster(b,Cluster.b);
                corpus[d].ViewCluster(c,Cluster.c);

                auto b_left_map= Cluster.b.first;

                for(auto it_b_left = b_left_map.begin(); it_b_left != b_left_map.end(); ++it_b_left){
                    Cluster.a_sx.insert(std::make_pair( it_b_left-> first , it_b_left-> second ));
                }

                auto c_left_map= Cluster.c.first;
                for(auto it_c_left = c_left_map.begin(); it_c_left != c_left_map.end(); ++it_c_left){

                    if(!(Cluster.a_sx.insert(std::make_pair( it_c_left-> first , it_c_left-> second )).second))
                        Cluster.a_sx[it_c_left-> first] = Cluster.b.first[it_c_left-> first] + it_c_left-> second;
                }

                auto b_right_map= Cluster.b.second;

                for(auto it_b_right = b_right_map.begin(); it_b_right != b_right_map.end(); ++it_b_right){
                    Cluster.a_dx.insert(std::make_pair( it_b_right-> first , it_b_right-> second ));
                }

                auto c_right_map= Cluster.c.second;

                for(auto it_c_right = c_right_map.begin(); it_c_right != c_right_map.end(); ++it_c_right){

                    if(!(Cluster.a_dx.insert(std::make_pair( it_c_right-> first , it_c_right-> second )).second))
                    Cluster.a_dx[it_c_right-> first] = Cluster.b.second[it_c_right-> first] + it_c_right-> second;
                }

                Cluster.a = std::make_pair(Cluster.a_sx,Cluster.a_dx);

                nj.b= corpus[d].ViewDataCountID(b);
                nj.c= corpus[d].ViewDataCountID(c);
                nj.a=nj.b+nj.c;

                nj.b_sub.first = corpus[d].ViewDataCountLeftID(b);
                nj.b_sub.second = corpus[d].ViewDataCountRightID(b);
                nj.c_sub.first = corpus[d].ViewDataCountLeftID(c);
                nj.c_sub.second = corpus[d].ViewDataCountRightID(c);
                nj.a_sub.first= nj.b_sub.first + nj.c_sub.first;
                nj.a_sub.second= nj.b_sub.second + nj.c_sub.second;

                corpus[d].InsertNewCluster(Cluster.a,pi.a,0.0,0.0,nj.a,nj.a_sub.first,nj.a_sub.second,0,0,0);

                Cluster.b.first.clear();
                Cluster.b.second.clear();
                Cluster.c.first.clear();
                Cluster.c.second.clear();
                Cluster.a_sx.clear();
                Cluster.a_dx.clear();
                Cluster.a.first.clear();
                Cluster.a.second.clear();
            }
        }


        wc.b.resize(D);
        wc.c.resize(D);
        wc.a.resize(D);

        unsigned int mb = 0;
        unsigned int mc = 0;

        #pragma omp parallel default(none) private(d,nj) shared(wc,b,c,mb,mc)
        {
            #pragma omp for schedule(static,1) reduction(+ : mb, mc)
            for(d=0; d<D;++d){
                nj.b= corpus[d].ViewDataCountID(b);
                nj.c= corpus[d].ViewDataCountID(c);

                mb += FindBestNumTable(Alpha, K, nj.b, LogStirlingNumbers);
                mc += FindBestNumTable(Alpha, K, nj.c, LogStirlingNumbers);
                wc.b[d]= nj.b;
                wc.c[d]= nj.c;
                wc.a[d]= nj.b + nj.c;
		  }
        }


        m.Tilde_b_c[0]=mb;
        m.Tilde_b_c[1]=mc;

        LogLikelihood_cap = this -> Model.Marginalized_Loglikelihood(NewClusterID);
        LogLikelihood_bc = this -> Model.Marginalized_Loglikelihood(b) + this -> Model.Marginalized_Loglikelihood(c);

        LogQ = static_cast<long double>(logL[b*K+b]) + static_cast<long double>(logL[c*K+c]) - static_cast<long double>(logL[b*K+c]) - static_cast<long double>(logL[c*K+b]);

        LogH += lgamma(static_cast<long double>( m.Tilde_b_c[0] + m.Tilde_b_c[1]))
                - std::log(static_cast<long double>(Gamma)) - lgamma(static_cast<long double>( m.Tilde_b_c[0])) - lgamma(static_cast<long double>( m.Tilde_b_c[1]));

        LogH += static_cast<long double>(m.Tilde_b_c[0])*std::log(static_cast<long double>(Beta.b_c[0])) + static_cast<long double>(m.Tilde_b_c[1])*std::log(static_cast<long double>(Beta.b_c[1]))
                - static_cast<long double>(m.Tilde_b_c[0] + m.Tilde_b_c[1])* std::log(static_cast<long double>(Beta.a));

        LogH += static_cast<long double>(LogLikelihood_cap) - static_cast<long double>(LogLikelihood_bc);

        LogH += LogQ;

        LogH += static_cast<long double>(logQ_split) - static_cast<long double>(logQ_merge);

            for(d=0; d<D; ++d){
                nj.a = wc.a[d];
                nj.b = wc.b[d];
                nj.c = wc.c[d];

                if (nj.a>0) LogH += lgamma(static_cast<long double>(Alpha*Beta.a + nj.a)) - lgamma(static_cast<long double>(Alpha*Beta.a));
                if (nj.b>0) LogH += lgamma(static_cast<long double>(Alpha*Beta.b_c[0])) - lgamma(static_cast<long double>(Alpha*Beta.b_c[0] + nj.b));
                if (nj.c>0) LogH += lgamma(static_cast<long double>(Alpha*Beta.b_c[1])) - lgamma(static_cast<long double>(Alpha*Beta.b_c[1] + nj.c));
            }


        LogH = std::min(LogH , static_cast<long double>(0));

        Gen.setNumThreads(1);

        u = Gen.runif();



        if(log(static_cast<long double>(u)) < LogH ){


            ProposedClusters.push_back(NewClusterID);
            this-> Gibbs_SubCluster (ProposedClusters);

            ClustersToBeRemoved.push_back(b);
            ClustersToBeRemoved.push_back(c);
            ++ NewClusterID;

		}
		else{
            Model.RemoveOneCluster(NewClusterID);

            omp_set_num_threads(OMP_NUM_THREADS);

            #pragma omp parallel for schedule (static,1) shared(NewClusterID)
            for(typename Corpus::size_type d=0; d<D;++d)
                corpus[d].RemoveCluster(NewClusterID);

        }

        wc.b.clear();
        wc.c.clear();
        wc.a.clear();
        LogH = 0.0;
        ProposedClusters.clear();


    }


    if(!(ClustersToBeRemoved.empty()))
		 Model.RemoveClusters(ClustersToBeRemoved);

    K = Model.ViewK();


    if(!(ClustersToBeRemoved.empty())){

	   omp_set_num_threads(OMP_NUM_THREADS);

	   #pragma omp parallel for schedule(static,1) shared(ClustersToBeRemoved)
	   for(typename Corpus::size_type d=0; d<D; ++d){
	      corpus[d].RemoveCluster(ClustersToBeRemoved);
           }

	}

}

//Calcolo di matrice 8.12, vedi Relazione_Parisi_Perego.pdf
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT, unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::computeLogL(){


    logL.reserve(K*K);
    for(size_t i=0; i<K*K; ++i)
        logL.push_back(0.0);

	vector<pair<POINT,unsigned int>> nidjk;
	double pijh = 0.0;
	double pijl = 0.0;

	omp_set_num_threads(OMP_NUM_THREADS);

	size_t h;
	size_t l;


	#pragma omp parallel for schedule(static,1) default(none) private(pijh,pijl,h,l,nidjk)
	for(h=0; h<K; ++h){
		for(l=0; l<K; ++l){

			if(h==l){
			    for(typename Corpus::size_type d=0; d<D;++d){
					corpus[d].ViewIdCounts(nidjk,h);
					pijh = corpus[d].ViewPiID(h);
					for(auto it = nidjk.cbegin(); it!= nidjk.cend(); ++it)
						logL[h*K + l] += (it->second)*std::log(pijh) + (it->second)*Model.Loglikelihood((it->first),h);

					nidjk.clear();
			   }
			} else{
			    for(typename Corpus::size_type d=0; d<D;++d){
			        corpus[d].ViewIdCounts(nidjk,h);
				    pijh = corpus[d].ViewPiID(h);
					pijl = corpus[d].ViewPiID(l);
					for(auto it = nidjk.cbegin(); it!= nidjk.cend(); ++it)
					    logL[h*K + l] +=  (it->second)* std::log( pijh*std::exp(Model.Loglikelihood((it->first),h)) + pijl*std::exp(Model.Loglikelihood((it->first),l)) );

					nidjk.clear();
				}

		    }
	    }
    }


}

//Calcola la quantità 3.25 della Relazione_Parisi_Perego.pdfCalcola la quantità 3.25 della Relazione_Parisi_Perego.pdf
template <template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT, unsigned int DIM>
long double HDP_MCMC<MODEL,DOCUMENT,DIM>::logq(){

    long double logq = 0.0;
    vector<pair<POINT,unsigned int>> nidjk;

    omp_set_num_threads(OMP_NUM_THREADS);
    double num;
    double den;

    #pragma omp parallel for private(num,den,nidjk) schedule(dynamic,1) reduction(+:logq)
    for(size_t k=0; k<K; ++k){
        for(typename Corpus::size_type d=0; d<D;++d){
            corpus[d].ViewIdCounts(nidjk,k);
            num = 0.0;

		    for(auto i: nidjk){
                den = 0.0;
				if(i.second>0){

					for(size_t l=0; l<K; ++l){
						den += corpus[d].ViewPiID(l) * std::exp( Model.Loglikelihood(i.first,l));
					}

					num = std::log(corpus[d].ViewPiID(k)) +  Model.Loglikelihood(i.first,k);
					logq += i.second *(num - std::log(den));
				}
            }

			nidjk.clear();
        }
	}

    return logq;

}

//Split Locale
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::LocalSplit(){

    CLUSTER Temp;
    POINT curr;
    BETA _Beta;

    NUMTABLE mjk;
    mjk.Tilde_b_c.assign(2,0);
    unsigned int NumTable_b;
    unsigned int NumTable_c;

    NJK njk;
    THETA _Theta;
    C _c;
    PI _Pi;

    vector<double> params;
    ClusterID Kcurrent;

	vector<ClusterID> NewClusterID;
    NewClusterID.resize(2);

    vector<ClusterID> K_temp;
    long double LogH;;
    long double Log_q_zb_cap;
    long double Log_q_zc_cap;
    double MargLogLike_a;
    double MargLogLike_b;
    double MargLogLike_c;

    bool DoSplit;

    double Q = (2 * floor(static_cast<double>(K + 1)/2))/(K*(K+1));

    double u;

    for (size_t k=0; k<K; k++){
        Kcurrent = Model.ViewK();
        DoSplit = true;

        DoSplit = this->IsEmptySubcluster(k);

        if(!DoSplit){


            NewClusterID[0] = Kcurrent;
            NewClusterID[1] = Kcurrent+1;

            Model.AddOneCluster(NewClusterID[0]);
			Model.AddOneCluster(NewClusterID[1]);

            _c.a.clear();
            Model[k].ViewStatisticsLeft(_c.a);
            Model[NewClusterID[0]].SetStatistics(_c.a);
            Model[NewClusterID[0]].SetStatisticsLeft(_c.a);

            _c.a.clear();
            _c.a.assign(W,0);

            Model[NewClusterID[0]].SetStatisticsRight(_c.a);


            Model[k].ViewThetaLeft(_Theta);
            Model[NewClusterID[0]].SetTheta(_Theta);

            _c.a.clear();

            Model[k].ViewStatisticsRight(_c.a);
            Model[NewClusterID[1]].SetStatistics(_c.a);
            Model[NewClusterID[1]].SetStatisticsLeft(_c.a);

            _c.a.clear();
            _c.a.assign(W,0);

            Model[NewClusterID[1]].SetStatisticsRight(_c.a);

            Model[k].ViewThetaRight(_Theta);
            Model[NewClusterID[1]].SetTheta(_Theta);

            NumTable_b = 0;
            NumTable_c = 0;
            mjk.Tilde_b_c.clear();
            mjk.Tilde_b_c.resize(2,0);

            Log_q_zb_cap = 0.0;
            Log_q_zc_cap = 0.0;

            omp_set_num_threads(OMP_NUM_THREADS);
            #pragma omp parallel for schedule(static,1) \
            default(none) private(Temp,_Pi,curr,njk) shared(k,NewClusterID) \
            reduction(+: Log_q_zb_cap)reduction(+: Log_q_zc_cap) reduction(+: NumTable_b)reduction(+: NumTable_c)
            for(typename Corpus::size_type d = 0; d < D; d++){
                Temp.b.first.clear();
                Temp.c.first.clear();
                Temp.b.second.clear();
                Temp.c.second.clear();
                Temp.a.first.clear();
                Temp.a.second.clear();

                corpus[d].ViewCluster(k,Temp.a);

                Temp.b.first = Temp.a.first;
                Temp.c.first = Temp.a.second;
                Temp.b.second = Temp.a.first;
                Temp.c.second = Temp.a.second;

                for(auto it = Temp.b.second.begin(); it != Temp.b.second.end(); it++)
                    (*it).second = 0;

                for(auto it = Temp.c.second.begin(); it != Temp.c.second.end(); it++)
                    (*it).second = 0;

                _Pi.Tilde_b_c.clear();
                _Pi.Tilde_b_c.resize(2);

                _Pi.a = corpus[d].ViewPiID(k);

                _Pi.Tilde_b_c[0] = corpus[d].ViewPiLeftID(k);
                _Pi.Tilde_b_c[1] = corpus[d].ViewPiRightID(k);

                _Pi.Tilde_b_c[0] = _Pi.Tilde_b_c[0] * _Pi.a;
                _Pi.Tilde_b_c[1] = _Pi.Tilde_b_c[1] * _Pi.a;

               for(auto it = Temp.b.first.begin(); it != Temp.b.first.end(); it++){
                    curr = it->first;

                    if(Temp.b.first[curr]>0)
                        Log_q_zb_cap += (static_cast<long double>(Temp.b.first[curr])) *
                                        (log(static_cast<long double>(_Pi.Tilde_b_c[0])) +
                                        static_cast<long double>(Model.Loglikelihood(curr,NewClusterID[0])) -
                                        log(static_cast<long double>(_Pi.Tilde_b_c[0] * std::exp(Model.Loglikelihood(curr,NewClusterID[0])) + _Pi.Tilde_b_c[1] * std::exp(Model.Loglikelihood(curr,NewClusterID[1])))));

                    if(Temp.c.first[curr] >0)
                        Log_q_zc_cap += (static_cast<long double>(Temp.c.first[curr])) * (log(static_cast<long double>(_Pi.Tilde_b_c[1])) +
                                        static_cast<long double>(Model.Loglikelihood(curr,NewClusterID[1])) -
                                        log(static_cast<long double>(_Pi.Tilde_b_c[1] * std::exp(Model.Loglikelihood(curr,NewClusterID[1])) + _Pi.Tilde_b_c[0] * std::exp(Model.Loglikelihood(curr,NewClusterID[0])))));
                }

                njk.b = corpus[d].ViewDataCountLeftID(k); //NOMi
                njk.c = corpus[d].ViewDataCountRightID(k);//NOMI

                NumTable_b += FindBestNumTable(Alpha, K, njk.b, LogStirlingNumbers);
                NumTable_c += FindBestNumTable(Alpha, K, njk.c, LogStirlingNumbers);

                corpus[d].InsertNewCluster (Temp.b,0.0,0.0,0.0,njk.b,njk.b,0,0,0,0);
                corpus[d].InsertNewCluster (Temp.c,0.0,0.0,0.0, njk.c,njk.c,0,0,0,0);
            }

            mjk.Tilde_b_c[0] = NumTable_b;
            mjk.Tilde_b_c[1] = NumTable_c;

            Model[NewClusterID[0]].SetGlobalTable( mjk.Tilde_b_c[0]);
            Model[NewClusterID[1]].SetGlobalTable( mjk.Tilde_b_c[1]);

            params.clear();
            params.resize(2);
            params[0] = static_cast<double>(mjk.Tilde_b_c[0]);
            params[1] = static_cast<double>(mjk.Tilde_b_c[1]);

             _Beta.a = Model[k].ViewBeta();
            _Beta.b_c.clear();
            _Beta.b_c.resize(2);

            Gen.setNumThreads(1);
            Gen.rdirichlet(params,_Beta.b_c);

            _Beta.b_c[0] = _Beta.b_c[0] * _Beta.a;
            _Beta.b_c[1] = _Beta.b_c[1] * _Beta.a;

            Model[NewClusterID[0]].SetBeta(_Beta.b_c[0]);
            Model[NewClusterID[1]].SetBeta(_Beta.b_c[1]);

            Gen.setNumThreads(1);

            Model.UpdateOneThetaCluster(NewClusterID[0], Gen);
            Model.UpdateOneThetaCluster(NewClusterID[1], Gen);

            LogH = 0.0;

            MargLogLike_b = Model.Marginalized_Loglikelihood(NewClusterID[0]);
            MargLogLike_c = Model.Marginalized_Loglikelihood(NewClusterID[1]);
            MargLogLike_a = Model.Marginalized_Loglikelihood(k);

            LogH += log(static_cast<long double>(Gamma)) + lgamma(static_cast<long double>(mjk.Tilde_b_c[0])) + lgamma(static_cast<long double>(mjk.Tilde_b_c[1])) - lgamma(static_cast<long double>(mjk.Tilde_b_c[0] + mjk.Tilde_b_c[1]));

            LogH += static_cast<long double>(mjk.Tilde_b_c[0] + mjk.Tilde_b_c[1])*log(static_cast<long double>(_Beta.a))- static_cast<long double>(mjk.Tilde_b_c[0])*log(static_cast<long double>(_Beta.b_c[0])) - static_cast<long double>(mjk.Tilde_b_c[1])*log(static_cast<long double>(_Beta.b_c[1]));

            LogH += static_cast<long double>(MargLogLike_b) + static_cast<long double>(MargLogLike_c) - static_cast<long double>(MargLogLike_a);

            LogH -= (Log_q_zb_cap + Log_q_zc_cap);

            LogH += log(static_cast<long double>(Q));

            omp_set_num_threads(OMP_NUM_THREADS);

            #pragma omp parallel for default(none) private(njk) shared(_Beta,k) reduction(+:LogH) schedule(dynamic,1)
            for(typename Corpus::size_type d = 0; d < D; d++){

                njk.b = corpus[d].ViewDataCountLeftID(k);
                njk.c = corpus[d].ViewDataCountRightID(k);

                if(njk.b > 0 || njk.c > 0)  LogH += lgamma(static_cast<long double>(Alpha * _Beta.a))  - lgamma(static_cast<long double> (Alpha*_Beta.a + njk.b + njk.c));
                if(njk.b > 0)   LogH += lgamma(static_cast<long double>(Alpha * _Beta.b_c[0] + njk.b)) - lgamma(static_cast<long double>(Alpha * _Beta.b_c[0]));
                if(njk.c > 0)   LogH += lgamma(static_cast<long double>(Alpha * _Beta.b_c[1] + njk.c)) - lgamma(static_cast<long double>(Alpha * _Beta.b_c[1]));
            }

            LogH = std::min(LogH ,static_cast<long double>(0.0));

            Gen.setNumThreads(1);

            u = Gen.runif();

            if(log(static_cast<long double>(u)) < LogH ){

                K_temp.push_back(k);


                this->Gibbs_SubCluster(NewClusterID);
            }
            else{
            Model.RemoveClusters(NewClusterID);

            omp_set_num_threads(OMP_NUM_THREADS);
            #pragma omp parallel for schedule (static,1)
            for(typename Corpus::size_type d=0; d<D; d++){
                corpus[d].RemoveCluster(NewClusterID);
            }
        }
        }
    }

    if(K_temp.size() != 0){

        omp_set_num_threads(OMP_NUM_THREADS);

        #pragma omp parallel for schedule(static,1)
       	for(typename Corpus::size_type d =0; d<D; d++){
            corpus[d].RemoveCluster(K_temp);
        }

        Model.RemoveClusters(K_temp);

    	K = Model.ViewK();

        this->UpdateDocWeights();
    }

    K = Model.ViewK();

}

//Global Split
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::GlobalSplit(){

    ClusterID K_proposed;

    vector<ClusterID> NewClusterID;
    NewClusterID.resize(2);
    NewClusterID[0] = K;
    NewClusterID[1]= K + 1;

    POINT curr;
    vector<ClusterID> K_current;
    bool DoSplit = true;
    unsigned int Index;

    Corpus OldCorpus;
    MODEL<DIM> OldModel;

    NUMTABLE mjk;
    unsigned int Thread_mjka;

    NJK njk;
    DATACOUNT DataCount;
    CLUSTER Cluster;
    THETA Theta;
    C _c;

    PI Pi;
    vector<double> Pi_a;

    vector<long double> ThreadLogH(OMP_NUM_THREADS,0);
    long double LogH = 0.0;
    long double LogQ = log(static_cast<long double>(2.0)) - log(static_cast<long double>(K+1));

    double u;

    OldCorpus = corpus;
    OldModel = Model;

	Model.ViewKey(K_current);

	Gen.setNumThreads(1);

    while(K_current.size()!=0){
        Index = Gen.runifdiscrete(K_current.size());
        DoSplit = this->IsEmptySubcluster(K_current[Index]);

        if(!DoSplit){
            K_proposed = K_current[Index];
            break;
        }
        else{
            K_current.erase(K_current.begin() + Index);
        }
    }


    if(DoSplit)return;


    vector<unsigned int> Thread_var;
    Thread_var.resize(K);

    omp_set_num_threads(OMP_NUM_THREADS);
    Thread_mjka = 0;

    #pragma omp parallel for private(mjk) default(none) shared(Thread_var) reduction(+:Thread_mjka) schedule(static,1)
    for(size_t k=0; k<K; k++){

        mjk.a = 0;
        for(typename Corpus::size_type d = 0; d<D; d++){
            mjk.a += FindBestNumTable(Alpha, K,  corpus[d].ViewDataCountID(k), LogStirlingNumbers);
        }

        Thread_var[k] = mjk.a;
        Thread_mjka += mjk.a;
    }

    mjk.Tilde_k = Thread_var;
    mjk.a = Thread_mjka;

    LogH += log(static_cast<long double>(Gamma)) + lgamma(static_cast<long double>(Gamma) + static_cast<long double>(mjk.a));

	for(size_t k =0; k<K; k++)
  		LogH -= static_cast<long double>(Model.Marginalized_Loglikelihood(k));

    #pragma omp parallel for default(none) shared(mjk) reduction(+ : LogH) schedule(static,1)
    for(size_t k=0; k<K; k++){
        LogH += static_cast<long double>(mjk.Tilde_k[k]) * log(static_cast<long double>(Model[k].ViewBeta())) - lgamma(static_cast<long double>(mjk.Tilde_k[k]));
    }

	#pragma omp parallel for default(none) schedule(dynamic,1) private(DataCount) reduction(+:LogH)
    for(typename Corpus::size_type d=0; d<D; d++){
        corpus[d].ViewDataCount(DataCount.a);

        for(size_t k=0; k<K; k++)
            if(DataCount.a[k] > 0) LogH += lgamma(static_cast<long double>(Alpha *Model[k].ViewBeta())) - lgamma(static_cast<long double>(Alpha *Model[k].ViewBeta() + DataCount.a[k]));

	}


    LogH += Model.LogDensity(K_proposed);

    LogH += logq();

	Model.AddOneCluster(NewClusterID[0]);
	Model.AddOneCluster(NewClusterID[1]);

    Model[K_proposed].ViewStatisticsLeft(_c.a);
    Model[NewClusterID[0]].SetStatistics(_c.a);
    Model[NewClusterID[0]].SetStatisticsLeft(_c.a);

    _c.a.clear();
    _c.a.assign(W,0);

    Model[NewClusterID[0]].SetStatisticsRight(_c.a);

    Model[K_proposed].ViewThetaLeft(Theta);
    Model[NewClusterID[0]].SetTheta(Theta);

    Model[K_proposed].ViewStatisticsRight(_c.a);
    Model[NewClusterID[1]].SetStatistics(_c.a);
    Model[NewClusterID[1]].SetStatisticsLeft(_c.a);

    _c.a.clear();
    _c.a.assign(W,0);

    Model[NewClusterID[1]].SetStatisticsRight(_c.a);

    Model[K_proposed].ViewThetaRight(Theta);
    Model[NewClusterID[1]].SetTheta(Theta);

    Model.RemoveOneCluster(K_proposed);

    #pragma omp parallel for default(none) shared(K_proposed) private(Cluster, Pi, njk) schedule(static,1)
	for(typename Corpus::size_type d=0; d<D; d++){

        Cluster.b.first.clear();
        Cluster.b.second.clear();
        Cluster.c.first.clear();
        Cluster.c.second.clear();

        corpus[d].ViewCluster(K_proposed, Cluster.a);

        Cluster.b.first = Cluster.a.first;
        Cluster.c.first = Cluster.a.second;
        Cluster.b.second = Cluster.a.first;
        Cluster.c.second = Cluster.a.second;

        for(auto it = Cluster.b.second.begin(); it != Cluster.b.second.end(); it++)
            (*it).second = 0;

        for(auto it = Cluster.c.second.begin(); it != Cluster.c.second.end(); it++)
            (*it).second = 0;

        Pi.a = corpus[d].ViewPiID(K_proposed);

        Pi.b_c.clear();
        Pi.b_c.resize(2);
        Pi.b_c[0] = Pi.a * corpus[d].ViewPiLeftID(K_proposed);
        Pi.b_c[1] = Pi.a * corpus[d].ViewPiRightID(K_proposed);

        njk.b = corpus[d].ViewDataCountLeftID(K_proposed);
        njk.c = corpus[d].ViewDataCountRightID(K_proposed);

        corpus[d].InsertNewCluster(Cluster.b,Pi.b_c[0],0.0,0.0, njk.b,njk.b,0,0,0,0);
        corpus[d].InsertNewCluster(Cluster.c,Pi.b_c[1],0.0,0.0, njk.c,njk.c,0,0,0,0);

        corpus[d].RemoveCluster(K_proposed);

    }

    K = Model.ViewK();

    this->UpdateAssignment_Cluster();

    LogH -= this->logq();

    #pragma omp parallel for default(none) private(Thread_mjka) schedule(static,1)
    for(size_t k=0; k<K; k++){
		Thread_mjka = 0;

        for(typename Corpus::size_type d=0; d<D; d++){
            Thread_mjka += FindBestNumTable(Alpha, K, corpus[d].ViewDataCountID(k), LogStirlingNumbers);
        }

       Model[k].SetGlobalTable(Thread_mjka);

    }

    this->UpdateBeta();

    this->UpdateDocWeights();

    Model.UpdateThetaCluster(Gen);

    //omp_set_num_threads(OMP_NUM_THREADS);

    Thread_mjka = 0;

    #pragma omp parallel for reduction(+: Thread_mjka) schedule(static,1)
    for(size_t k=0; k<K ;k++){
        Thread_mjka += Model[k].ViewGlobalTable();
    }

    LogH += -lgamma(static_cast<long double>(Gamma + Thread_mjka));

    for(size_t k=0; k<K; k++)
        LogH += static_cast<long double>(Model.Marginalized_Loglikelihood(k));

    LogH += LogQ;

    #pragma omp parallel for default(none) reduction(+: LogH) private(mjk) schedule (static,1)
    for(size_t k=0; k<K; k++){
        mjk.a = Model[k].ViewGlobalTable();
        LogH += lgamma(static_cast<long double>(mjk.a)) - static_cast<long double>(mjk.a)*(Model[k].ViewBeta());
    }

    #pragma omp parallel for default(none) schedule (dynamic,1) reduction(+:LogH) private(DataCount)
    for(typename Corpus::size_type d=0; d<D; d++){
        DataCount.a.clear();
        corpus[d].ViewDataCount(DataCount.a);

        for(size_t k=0; k<K; k++)
            if(DataCount.a[k]>0)LogH += lgamma(static_cast<long double>(Alpha) * Model[k].ViewBeta() + static_cast<long double>(DataCount.a[k])) - lgamma(static_cast<long double>(Alpha) * Model[k].ViewBeta());
    }

    LogH = std::min(LogH ,static_cast<long double>(0.0));

    Gen.setNumThreads(1);

	u = Gen.runif();

    if(log(u) < LogH ){

        NewClusterID.clear();

        Model.ViewKey(NewClusterID);

        this->Gibbs_SubCluster(NewClusterID);
    }
	else{
        this->Swap(OldCorpus,corpus);

        this->Swap(OldModel,Model);

        K = Model.ViewK();
    }


}

//Global Merge
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::GlobalMerge(){

    if(K==1)
        return;

    vector<ClusterID> CurrentClusters;
    Model.ViewKey(CurrentClusters);

	Gen.setNumThreads(1);

	unsigned int guess = Gen.runifdiscrete(K);
	ClusterID b = CurrentClusters[guess];
	CurrentClusters.erase(CurrentClusters.begin() + guess);
	guess = Gen.runifdiscrete(K-1);
	ClusterID c = CurrentClusters[guess];

	ClusterID NewClusterID = K;

	C stats;
	NJK nj;
	NUMTABLE m;
	PI pi;
	CLUSTER Cluster;
	BETA Beta;

	Corpus Oldcorpus = corpus;
	MODEL<DIM> OldModel= Model;

	long double logq = 0.0;
	double Loglikelihood = 0.0;
	long double LogDensity = 0.0;

	double logQsplit = -log(K-1);
	double logQmerge = log(2) - log(K) - log(K-1);

	long double LogH = 0.0;

    double u = 0.0;

	Model.AddOneCluster(NewClusterID);

	Model[b].ViewStatistics(stats.b);
	Model[c].ViewStatistics(stats.c);

	if(stats.b.size()!= stats.c.size()){
	  std::cerr<< "In Global Merge, c.b and c.c must have the same dimension"<<std::endl;
	  exit(1);
	}


	Model[NewClusterID].UpdateStatistics(stats.b,stats.c);

	Gen.setNumThreads(1);

    Model.UpdateOneThetaCluster(NewClusterID,Gen);

    m.Tilde_k.resize(K,0);
	unsigned int mtildesum = 0;
	unsigned int mtildek = 0;
	size_t k;
	typename Corpus::size_type d;

	omp_set_num_threads(OMP_NUM_THREADS);

    #pragma omp parallel default(none) shared(mtildesum,m) private(nj,k,d) firstprivate(mtildek)
	{
        #pragma omp for schedule(static,1) reduction(+ : mtildesum)
        for(k=0; k<K; ++k){
            for(d=0;d<D;++d){
                nj.k = corpus[d].ViewDataCountID(k);//NM
                mtildek += FindBestNumTable(Alpha, K, nj.k, LogStirlingNumbers);
            }


            m.Tilde_k[k] = mtildek;
            mtildesum += mtildek;
            mtildek = 0;
        }
	}

	m.Tilde_sum = mtildesum;

	logq = this-> logq();

	for(size_t k=0; k<K; ++k)
	   Loglikelihood += Model.Marginalized_Loglikelihood(k);

	LogH += lgamma(static_cast<long double>( Gamma + m.Tilde_sum)) - log(static_cast<long double>(Gamma));

	LogH += logq - static_cast<long double>(Loglikelihood);

	LogH +=  static_cast<long double>(logQsplit) - static_cast<long double>(logQmerge);


    #pragma omp parallel for default(none) private(k,Beta) shared(m) schedule(static,1) reduction (+:LogH)

		for(k=0; k<K; ++k){
            Beta.k = Model[k].ViewBeta();
            LogH += static_cast<long double>(m.Tilde_k[k])*log(static_cast<long double>(Beta.k)) - lgamma(static_cast<long double>(m.Tilde_k[k]));
	    }

    #pragma omp parallel for default(none) private(d,k,Beta,nj) schedule(dynamic,1) reduction(+:LogH)

		for(d=0;d<D;++d){
			for(k=0;k<K;++k){
			    Beta.k = Model[k].ViewBeta();
				nj.k = corpus[d].ViewDataCountID(k);//nome
				if(nj.k > 0) LogH+= lgamma(static_cast<long double>(Alpha*Beta.k))-lgamma(static_cast<long double>(Alpha*Beta.k + nj.k));
			}
		}

    LogDensity = Model.LogDensity(NewClusterID);

	Model.RemoveOneCluster(b);
	Model.RemoveOneCluster(c);

	K = Model.ViewK();

	pi.b_c.resize(2);
	pi.a = 0.0;

	#pragma omp parallel default(none) private(d) shared(b,c,Cluster) firstprivate(pi)
	{
        #pragma omp for schedule(static,1)
        for(d =0; d<D; ++d){
            pi.b_c[0] = corpus[d].ViewPiID(b); //pib
            pi.b_c[1] = corpus[d].ViewPiID(c); //pic
            pi.a = pi.b_c[0] + pi.b_c[1];

            corpus[d].InsertNewCluster(Cluster.a,pi.a,0.0,0.0,0,0,0,0,0,0);

            corpus[d].RemoveCluster(b,c);
        }
	}

	this-> UpdateAssignment_Cluster();

	logq = this-> logq();

	mtildesum = 0;
	mtildek = 0;

	for(auto& i: m.Tilde_k )
	    i = 0;


    #pragma omp parallel default(none) shared(mtildesum,m) private(nj,k,d) firstprivate(mtildek)
	{
	    #pragma omp for schedule(static,1) reduction(+ : mtildesum)
        for(k=0; k<K; ++k){

                for(d=0;d<D;++d){
                    nj.k = corpus[d].ViewDataCountID(k);//nome
                    mtildek += FindBestNumTable(Alpha, K, nj.k, LogStirlingNumbers);
                }


            m.Tilde_k[k] = mtildek;
            mtildesum += mtildek;
            mtildek = 0;

            Model[k].SetGlobalTable(m.Tilde_k[k]);
        }
    }

	m.Tilde_sum = mtildesum;

	this -> UpdateBeta();

	this -> UpdateDocWeights();

	Model.UpdateThetaCluster(Gen);

    LogH -= lgamma(static_cast<long double>( Gamma + m.Tilde_sum));

    Loglikelihood = 0.0;

	for(size_t k=0; k<K; ++k)
	   Loglikelihood += Model.Marginalized_Loglikelihood(k);

	LogH += static_cast<long double>(Loglikelihood) - logq - LogDensity;

    #pragma omp parallel for default(none) private(k,Beta) shared(m) schedule(static,1) reduction (+:LogH)
	for(k=0; k<K; ++k){
        Beta.k = Model[k].ViewBeta();
        LogH += lgamma(static_cast<long double>(m.Tilde_k[k])) - static_cast<long double>(m.Tilde_k[k])*log(static_cast<long double>(Beta.k));
	}

    #pragma omp parallel for default(none) private(d,k,Beta,nj) schedule(dynamic,1) reduction(+:LogH)
	for(d=0;d<D;++d){
		for(k=0;k<K;++k){
			Beta.k = Model[k].ViewBeta();
			nj.k = corpus[d].ViewDataCountID(k); //nome
			if(nj.k > 0) LogH+= lgamma(static_cast<long double>(Alpha*Beta.k + nj.k))- lgamma(static_cast<long double>(Alpha*Beta.k));
		}
	}

	LogH = std::min(LogH , static_cast<long double>(0));

	Gen.setNumThreads(1);

	u = Gen.runif();

	if(log(static_cast<long double>(u)) < LogH ){

        vector<ClusterID> ProposedClusters;

        Model.ViewKey(ProposedClusters);

        this-> Gibbs_SubCluster(ProposedClusters);
     }
     else{
        this -> Swap (Oldcorpus,corpus);

        this -> Swap(OldModel,Model);

        K = Model.ViewK();
    }

}

template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SaveRunTime(){

    std::ofstream fileK("./cpp_results/AllK1.txt", std::ios::app);

	if (fileK.fail()){
		std::cerr <<"ERROR GENERATED (SaveAllK)!" <<std::endl;
		std::cerr <<"Cannot open AllAlpha.txt" <<std::endl;
		exit(1);
	}

    std::ofstream fileAlpha("./cpp_results/AllAlpha1.txt", std::ios::app);

	if (fileAlpha.fail()){
		std::cerr <<"ERROR GENERATED (SaveAllAlpha)!" <<std::endl;
		std::cerr <<"Cannot open AllAlpha.txt" <<std::endl;
		exit(1);
	}

	std::ofstream fileGamma("./cpp_results/AllGamma1.txt", std::ios::app);

	if (fileGamma.fail()){
		std::cerr <<"ERROR GENERATED (SaveAllGamma)!" <<std::endl;
		std::cerr <<"Cannot open AllAlpha.txt" <<std::endl;
		exit(1);
	}

	fileK<<K<<std::endl;
	fileK.close();

	fileAlpha<<Alpha<<std::endl;
	fileAlpha.close();

	fileGamma<<Gamma<<std::endl;
	fileGamma.close();

}
//Salvataggio su file delllo storico dei Clusters
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SaveAllK(){

	std::ofstream file("./cpp_results/AllK.txt");

	if (file.fail()){
		std::cerr <<"ERROR GENERATED (SaveAllAlpha)!" <<std::endl;
		std::cerr <<"Cannot open AllAlpha.txt" <<std::endl;
		exit(1);
	}

	for (vector<unsigned int>::const_iterator it=AllK.cbegin(); it!=AllK.cend(); it++)
		file <<(*it) <<std::endl;

	file.close();

}

//Salvataggio su file dello storico degli Alpha
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SaveAllAlpha(){

	std::ofstream file("./cpp_results/AllAlpha.txt");

	if (file.fail()){
		std::cerr <<"ERROR GENERATED (SaveAllAlpha)!" <<std::endl;
		std::cerr <<"Cannot open AllAlpha.txt" <<std::endl;
		exit(1);
	}

	for (vector<double>::const_iterator it=AllAlpha.cbegin(); it!=AllAlpha.cend(); it++)
		file <<(*it) <<std::endl;

	file.close();
}

//Slavataggio su file dello storico dei Gamma
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SaveAllGamma(){

	std::ofstream file("./cpp_results/AllGamma.txt");

	if (file.fail()){
		std::cerr <<"ERROR GENERATED (SaveAllGamma)!" <<std::endl;
		std::cerr <<"Cannot open AllGamma.txt" <<std::endl;
		exit(1);
	}

	for (vector<double>::const_iterator it=AllGamma.cbegin(); it!=AllGamma.cend(); it++)
		file <<(*it) <<std::endl;

	file.close();
}

//Salvataggio su file dei pesi globali dei cluster delle ultime 100 iterazioni
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SaveLastBeta(const std::string& Filename){

    std::ofstream file (Filename, std::ios::binary| std::ios::app);

    if (file.fail()){
		std::cerr <<"ERROR GENERATED (SaveLast100Beta)!" <<std::endl;
		std::cerr <<"Cannot open LastBeta.bin" <<std::endl;
		exit(1);
	}

	/*
	if((MaxIt>=100) && ((MaxIt - It == 100))){
        file.write(reinterpret_cast<char*>(&It), sizeof(unsigned int));
        file.flush();
	}else if((MaxIt<100) && ((It == 0))){
        file.write(reinterpret_cast<char*>(&It), sizeof(unsigned int));
        file.flush();
	}
	*/

	vector<double> B(K,0);
	unsigned int dim;

	for(size_t k=0; k<K; k++)
        B[k] = Model[k].ViewBeta();

    dim = B.size();

    //file.write(reinterpret_cast<char*>(&dim), sizeof(unsigned int));
    //file.flush();
    file.write(reinterpret_cast<char*>(B.data()),dim*sizeof(double));
    file.flush();
    file.close();

}



//Salvataggio su file dei pesi dei cluster specifici dei documenti delle ultime 100 iterazioni
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SaveLastPi(const std::string& FileName){

    std::ofstream file (FileName,std::ios::binary| std::ios::app);

    if (file.fail()){
		std::cerr <<"ERROR GENERATED (SaveLast100Pi)!" <<std::endl;
		std::cerr <<"Cannot open LastPi.bin" <<std::endl;
		exit(1);
	}

	/*
	if((MaxIt>=100) && ((MaxIt - It == 100))){
        file.write(reinterpret_cast<char*>(&It), sizeof(unsigned int));
        file.flush();
	}else if((MaxIt<100) && ((It == 0))){
        file.write(reinterpret_cast<char*>(&It), sizeof(unsigned int));
        file.flush();
	}
	*/

	vector<double> pi;
	unsigned int dim = K;

	//file.write(reinterpret_cast<char*>(&dim), sizeof(unsigned int));
    //file.flush();


	for(size_t d=0; d<D; ++d){
		corpus[d].ViewPi(pi);
		file.write(reinterpret_cast<char*>(pi.data()),(dim)*sizeof(double));
        file.flush();

	}

    file.close();

}

//Slvataggio su file dello dei parametri latenti delle ultime 100 iterazioni
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SaveLastTheta(const std::string& FileName){

    //std::ofstream file (FileName,std::ios::binary | std::ios::app );

   /*
   if (file.fail()){
		std::cerr <<"ERROR GENERATED (LastTheta)!" <<std::endl;
		std::cerr <<"Cannot open LastTheta" <<std::endl;
		exit(1);
	}*/

	/*if((MaxIt>=100) && ((MaxIt - It == 100))){
        file.write(reinterpret_cast<char*>(&It), sizeof(unsigned int));
        file.flush();
	}else if((MaxIt<100) && ((It == 0))){
        file.write(reinterpret_cast<char*>(&It), sizeof(unsigned int));
        file.flush();
	}*/

    //file.write(reinterpret_cast<char*>(&K), sizeof(unsigned int));
    //file.close();
    Model.PrintTheta(FileName);

}

//Salvataggio su file delle quantità necessarie per il calcolo del LPML
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::LPML(){

    vector<double> DataLikelihood;
    DataLikelihood.reserve(N);
    vector<pair<POINT,ClusterID>> Data_Cluster;

    POINT x_ji;
    ClusterID k;

    for(typename Corpus::size_type d=0; d<D; d++){
        corpus[d].ViewLabel(Data_Cluster);
        for(auto ItData=Data_Cluster.begin(); ItData != Data_Cluster.end(); ItData++){
            x_ji = (*ItData).first;
            k = (*ItData).second;
            DataLikelihood.push_back(std::exp(Model.Loglikelihood(x_ji,k)));
        }

        Data_Cluster.clear();
    }

    for(size_t n=0; n<N; n++)
        CPO[n] += DataLikelihood[n];

    if(It == MaxIt -1){
        std::ofstream LPML ("./cpp_results/CPO.bin",std::ios::binary | std::ios::app);
        if (LPML.fail()){
            std::cerr <<"ERROR GENERATED (LPML)!" <<std::endl;
            std::cerr <<"Cannot open CPO.bin" <<std::endl;
            exit(1);
        }
        LPML.write(reinterpret_cast<char*>(CPO.data()),CPO.size()*sizeof(double));
        LPML.close();
    }

}

//Slavataggio su file delle etichette per ogni dato ad ogni iterazione
template<template <unsigned int> class MODEL, template <unsigned int> class DOCUMENT ,unsigned int DIM>
void HDP_MCMC<MODEL,DOCUMENT,DIM>::SaveLabels(const std::string& FileName){

	vector<pair<unsigned int,ClusterID>> Data;
	Data.reserve(N);

	for(typename Corpus::size_type d=0; d<D;++d)
		corpus[d].ViewLabel(Data);

	if(Data.size()!=N){
		std::cerr<<"Error in SaveLabels: Data must have length N"<<std::endl;
		exit(1);
	}

	std::ofstream file (FileName, std::ios::binary | std::ios::app);   // app perchè ad ogni iterazione devo aggiungere il clustering
	//std::ofstream file2 ("OrdineDati.txt",std::ios::app);

	if(file.fail()){
		std::cerr<<"Cannot open Labels.bin"<<std::endl;
		exit(1);
    }
/*
    if(file2.fail()){
        std::cerr<<"non posso aprire il file"<<std::endl;
        exit(1);
    }
*/
	ClusterID label;
	for(vector<pair<unsigned int,ClusterID>>::const_iterator it = Data.cbegin(); it != Data.cend() ; ++it){

  //      file2<<it->first<<std::endl;
		label = it->second;
		file.write(reinterpret_cast<char*>(&label), sizeof(unsigned int));
	}
	//file2.close();
	file.close();

}


#endif
