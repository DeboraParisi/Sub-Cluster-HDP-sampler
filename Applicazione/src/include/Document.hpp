#ifndef _DOCUMENT_HPP
#define _DOCUMENT_HPP

#include <unordered_map>
#include <utility>
#include <vector>
#include "Functions.hpp"
#include "Type.hpp"

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <tuple>
#include <fstream>

using std::unordered_map;
using std::pair;
using std::vector;

/*! \file Document.hpp
 *
 * \brief In questo file sono presenti le classi che gestiscono i documenti o, più in generale, i gruppi di dati. 
 *  La classe generica fornisce l'interfaccia comune, mentre le classi derivate e specializzate sono specifiche del modello.
 * \date Febbraio 2016
 */


/*! \brief Classe generica per i gruppi
 *
 *  Classe astratta dove tutti i metodi virtuali sono null.
 *  Contiene i metodi che devono essere obbligatoriamente definiti in tutte le classi derivate.
 *  
 * \date Febbraio 2016
 */


template<typename Type, unsigned int DIM>
class GenericDocument{

	public:
        /*! \brief Aggiorna i pesi dei cluster specifici del gruppo; si veda equazione (3.5) in Relazione_Parisi_Perego.pdf
		 *  \param _AllBeta - pesi globali dei cluster
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdatePi (const vector<double>& , omprng&  ) = 0;

		/*! \brief Aggiorna i pesi, specifici del gruppo, dei subcluster del cluster k; si veda equazione (3.9) in Relazione_Parisi_Perego.pdf
		 *  \param _BetaLeft - peso globale del subcluster sinistro del cluster k
		 *  \param _BetaRight - peso globale dei subcluster destro del cluster k
		 *  \param k - id del cluster
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdatePiSub (const double , const double , const unsigned int , omprng& ) = 0;
		
				
	   	/*! \brief Aggiorna i pesi, specifici del gruppo, di tutti i subcluster; ; si veda equazione (3.9) in Relazione_Parisi_Perego.pdf
		 *  \param _BetaLeft - pesi globali dei subcluster sinistri
		 *  \param _BetaRight - pesi globali dei subcluster destri 
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdateAllPiSub(const vector<double> , const vector<double> , omprng& ) = 0;
		
		/*! \brief Aggiorna i tavoli; si veda equazione (3.3) in Relazione_Parisi_Perego.pdf
		 *  \param _stirling - numeri di stirling
		 *  \param _Beta - pesi globali dei cluster
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdateLocalTable (const vector<long double>& , const vector<double>& , omprng& ) = 0;
		
		/*! \brief Aggiorna i tavoli dei subcluster del cluster k; si veda equazione (3.12) in Relazione_Parisi_Perego.pdf
		 *  \param _stirling - numeri di stirling
		 *  \param _BetaLeft - peso globale del subcluster sinistro del cluster k
		 *  \param _BetaRight - peso globale del subcluster destro del cluster k
		 *  \param k - id del cluster
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdateLocalTableSub_OneCluster(const vector<long double>& , const double _, const double , const unsigned int, omprng& ) = 0;
		
		/*! \brief Aggiorna i tavoli dei subcluster di tutti i cluster; si veda equazione (3.12) in Relazione_Parisi_Perego.pdf
		 *  \param _stirling - numeri di stirling
		 *  \param _BetaLeft - pesi globali dei subcluster sinistri
		 *  \param _BetaRight - pesi globali dei subcluster destri 
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdateAllLocalTableSub (const vector<long double>& , const vector<double>& , const vector<double>& , omprng& ) = 0;
		
		 /*! \brief Aggiorna l'etichetta per il cluster di un dato, campionata con il metodo Sampling; si veda equazione (3.7) in Relazione_Parisi_Perego.pdf
		 *  \param _ThetaId - vettore dei pesi del dato in tutti i cluster
		 *  \param _VettId - id del dato
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdateZeta (const typename Type::THETA& , const unsigned int , omprng& )= 0;

		/*! \brief Aggiorna l'etichetta per il cluster e per il subcluster di un dato, campionate con il metodo Sampling; 
		 *  si vedano equazioni (3.7) - (3.11) in Relazione_Parisi_Perego.pdf
		 *  \param _ThetaId - vettore dei pesi del dato in tutti i cluster
		 *  \param _ThetaIdLeft - vettore dei pesi del dato in tutti i subcluster sinistri
		 *  \param _ThetaIdRight - vettore dei pesi del dato in tutti i subcluster destri
		 *  \param _VettId - id del dato
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdateZeta_and_Sub(const typename Type::THETA&, const typename Type::THETA& ,const typename Type::THETA& , const unsigned int , omprng& )= 0;

		/*! \brief Distribuisce il dato nei subcluster del cluster k, dopo aver campionato l'etichetta del subcluster con il metodo Sampling;
		 *   si veda equazione (3.11) in Relazione_Parisi_Perego.pdf
		 *  \param _ThetaIdLeft - vettore dei pesi del dato id in tutti i subcluster sinistri
		 *  \param _ThetaIdRight - vettore dei pesi del dato id in tutti i subcluster destri
		 *  \param id - id del dato
		 *  \param nidjk - numero di volte che il dato id nel gruppo j è capitata nel cluster k
		 *  \param k - id del cluster
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		virtual void UpdateZetaSub(const typename Type::THETA& ,const typename Type::THETA&  ,const typename Type::Point , const unsigned int , const unsigned int , omprng& ) = 0;

		/*! \brief Metodo che serve nelle mosse di M-H per rimuovere il cluster k
		 *  \param _k - id del cluster
		 */
		virtual void UpdateZeta(const unsigned int ) = 0;

		/*! \brief Metodo che serve nelle mosse di M-H per rimuovere due cluster
		 *  \param _k1 - id del cluster
		 *  \param _k2 - id del cluster
		 */
		virtual void UpdateZeta(const unsigned int , const unsigned int ) = 0;

		/*! \brief Visualizza il numero di dati nel gruppo j
		 *  \return Numero di dati nel gruppo j
		 */
		virtual unsigned int ViewNj() const = 0; 

	     /*! Estrae gli identificativi dei dati nel gruppi j
		 *  \param _VettId - viene riempito con gli id dei dati nel gruppo j
		 */
		virtual void ViewData( vector<typename Type::Point>& ) const = 0;      

		/*! \brief Visualizza i conteggi necessari per aggiornare i parametri latenti dei subcluster
		 *  \param _k - id del cluster
		 *  \param _counts4cleft - conteggi per parametro latente del subcluster sinistro
		 *  \param _counts4cright - conteggi per parametro latente del subcluster destro
		 */
		virtual void ViewCounts4c(const unsigned int , typename Type::STAT& , typename Type::STAT& )= 0; 

		/*! \brief Estrae il numero di tavoli in uno specifico cluster
		 *  \param _k - id del cluster
		 *  \return Numero di tavoli nel ristorante j che servono il piatto k
		 */
		virtual unsigned int ViewNumTableID(const unsigned int) const = 0;

		/*! \brief Estrae il numero di tavoli del subcluster sinistro di uno specifico cluster
		 *  \param _k - id del cluster
		 *  \return Numero di tavoli nel ristorante j che servono il piatto k left
		 */
		virtual unsigned int ViewNumTableLeftID(const unsigned int ) const = 0;

		/*! \brief Estrae il numero di tavoli del subcluster destro di uno specifico cluster
		 *  \param _k - id del cluster
		 *  \return Numero di tavoli nel ristorante j che servono il piatto k right
		 */
		virtual unsigned int ViewNumTableRightID(const unsigned int) const = 0;
		
		/*! \brief Estrae il numero di dati nel cluster k
		 *  \param _k - id del cluster
		 *  \return Numero di dati nel cluster k
		 */
		virtual unsigned int ViewDataCountID(const unsigned int ) const = 0;
		
		/*! \brief Estrae il numero di dati nel subcluster sinistro del cluster k
		 *  \param _k - id del cluster
		 *  \return Numero di dati nel subcluster sinistro del cluster k
		 */
		virtual unsigned int ViewDataCountLeftID(const unsigned int ) const = 0;

		/*! \brief Estrae il numero di dati nel subcluster destro del cluster k
		 *  \param _k - id del cluster
		 *  \return Numero di dati nel subcluster destro del cluster k
		 */
		virtual unsigned int ViewDataCountRightID(const unsigned int ) const = 0;
		
		/*! \brief Azzera il conteggio dei dati nel cluster k
		 *  \param k - id del cluster
		 */
		virtual void ResetDataCountSub(const unsigned int) = 0;
		
		/*! \brief Estrae il vettore del numero di dati in ogni cluster
		 *  \param _WordCount - vettore del numero di dati in ogni cluster
		 */
		virtual void ViewDataCount(vector<unsigned int>& ) const = 0;
		
		/*! \brief Estrae il vettore del numero di dati in ogni subcluster sinistro
		 *  \param _WordCountLeft - vettore del numero di dati in ogni subcluster sinistro
		 */
		virtual void ViewDataCountLeft(vector<unsigned int>& ) const = 0;

		/*! \brief Estrae il vettore del numero di dati in ogni subcluster destro
		 *  \param _WordCountRight - vettore del numero di dati in ogni subcluster destro
		 */
		virtual void ViewDataCountRight(vector<unsigned int>& ) const = 0;
		
		/*! \brief Estrae identificativi e conteggi dei dati nel cluster k
		 *  \param _nidjk - struttura che contiene identificativi e conteggi dei dati nel cluster k
		 *  \param _k - id del cluster
		 */
		virtual void ViewIdCounts(vector<pair<typename Type::Point,unsigned int>>& ,const unsigned int ) = 0; 
		
		/*! \brief Estrae il cluster k
		 *  \param _k - id del cluster da estrarre
		 *  \param _Cluster - struttura in cui viene estratto il cluster
		 */
		virtual void ViewCluster(const unsigned int, pair<unordered_map<typename Type::Point,unsigned int>,unordered_map<typename Type::Point,unsigned int>>& ) = 0;  
		
		/*! \brief Estrae il peso, specifico del gruppo, del cluster k
		 *  \param _k - id del cluster
		 *  \return Peso specifico del gruppo del cluster k
		 */
		virtual double ViewPiID(const unsigned int )const = 0;
		
		/*! \brief Estrae il peso, specifico del gruppo, del subcluster sinistro del cluster k
		 *  \param _k - id del topic
		 *  \return Peso specifico del gruppo del subcluster sinistro del cluster k
		 */
		virtual double ViewPiLeftID(const unsigned int ) const = 0;

		/*! \brief Estrae il peso, specifico del gruppo, del subcluster destro del cluster k
		 *  \param _k - id del cluster
		 *  \return Peso specifico del documento del subcluster destro del cluster k
		 */
		virtual double ViewPiRightID(const unsigned int) const = 0;

		/*! \brief Estrae il vettore di pesi dei cluster specifici dei gruppi
		 *  \param _pi - vettore di pesi dei cluster specifici dei gruppi
		 */
		virtual void ViewPi( vector<double>& ) const = 0;

		/*! \brief Estrae il vettore di pesi dei subcluster sinistri specifici dei gruppi
		 *  \param _pi - vettore di pesi dei subcluster sinistri specifici dei gruppi
		 */
		virtual void ViewPiLeft( vector<double>& ) const = 0;
		
		/*! \brief Estrae il vettore di pesi dei subcluster destri specifici dei gruppi
		 *  \param _pi - vettore di pesi dei subcluster destri specifici dei gruppi
		 */
		virtual void ViewPiRight( vector<double>& ) const = 0;

		/*! \brief Imposta il parametro di concentrazione del processo di Dirichlet che governa il gruppo
		 *  \param _alpha - parametro di concentrazione del processo di Dirichlet che governa il gruppo
		 */
		virtual void SetAlpha (const double ) = 0;

		/*! \brief Imposta il numero di dati nel gruppo j
		 *  \param _Nj - numero di dati nel gruppo j
		 */
		virtual void SetNj (const unsigned int ) = 0;
		
		/*! \brief Imposta il vettore di pesi dei cluster specifici dei gruppi
		 */
		virtual void SetPi(vector<double>& ) = 0;

	    /*! \brief Inserisce un nuovo cluster
		 *  \param NewCluster - il nuovo cluster
		 *  \param _Pi - peso del nuovo cluster specifico del gruppo
		 *  \param _PiLeft - peso del subcluster sinistro nuovo cluster, specifico del gruppo
		 *  \param _PiRight - peso del subcluster destro del nuovo cluster, specifico del gruppo
		 *  \param _WordCount - numero di dati nel nuovo cluster
		 *  \param _WordCountLeft - numero di dati nel subcluster sinistro del nuovo cluster
		 *  \param _WordCountRight - numero di dati nel subcluster destro del nuovo cluster
		 *  \param _LocalTable - numero di tavoli che servono il nuovo piatto nel ristorante j
		 *  \param _LocalTableLeft - numero di tavoli che servono il nuovo piatto left nel ristorante j
		 *  \param _LocalTable - numero di tavoli che servono il nuovo piatto right nel ristorante j
		 */
		virtual void InsertNewCluster(const pair < unordered_map < typename Type::Point, unsigned int>, unordered_map < typename Type::Point, unsigned int> >& , const double, const double, const double ,
	                   					   const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int )= 0;

		/*! \brief Rimuove i cluster con identificativo presente nel vettore in ingresso
		 *  \param _k - vettore con gli id dei cluster da eliminare
		 */
		virtual void RemoveCluster(const vector<unsigned int>& ) = 0;
		
		/*! \brief Rimuove un cluster
		 *  \param _k - id del cluster da eliminare
		 */
		virtual void RemoveCluster( const unsigned int ) = 0;
	
		/*! \brief Rimuove due cluster
		 *  \param _k1 - id del cluster da eliminare
		 *  \param _k2 - id del cluster da eliminare
		 */
		virtual void RemoveCluster(unsigned int ,unsigned int ) = 0;

		/*! \brief Verifica se un cluster ha il subcluster sinistro vuoto
		 *  \param _k - id del cluster di cui controllare il subcluster
		 */
		virtual unsigned int CheckLeftSubcluster(const unsigned int ) = 0; 

		/*! \brief Verifica se un cluster ha il subcluster destro vuoto
		 *  \param _k - id del cluster di cui controllare il subcluster
		 */
		virtual unsigned int CheckRightSubcluster(const unsigned int ) = 0;  

		/*! \brief Estrae le etichette associate ai dati
		 *  \param Data - struttura in cui estrarre le etichette
		 */
		virtual void ViewLabel(vector<pair<typename Type::Point,unsigned int>>& ) = 0;  

		//virtual void PrintData() = 0;  
	
		/*! \brief Acquisisce i dati
		 *  \param SSTR - contiene id del dato e numero di volte che il dato compare nel gruppo
		 */
		virtual void SetDataset(std::istringstream&) = 0;
	
        /*! \brief Smista i dati nel contenitore Zeta
		 *  \param _K - numero iniziale di cluster
		 *  \param Gen - generatore di numeri casuali
		 */	
		virtual unsigned int SortData(unsigned int, omprng&) = 0;


    private:

	   	/*! \brief Aggiorna i conteggi dei dati nei cluster
		 */
		virtual void UpdateDataCount() = 0;
	
		/*! \brief Campionamento da distribuzione categorica per l'etichetta del cluster o del subcluster
		 *  \param _temp_counts - vettore che contiene i conteggi del dato id nei cluster
		 *  \param _Weights - pesi con cui campionare le/la etichette/a
		 *  \param _nidj - numero di volte che il dato id compare nel gruppo j
		 *  \param Gen - generatore di numeri casuali
		 */
		virtual void Sampling (std::vector<unsigned int>& , std::vector<double>& , unsigned int, omprng& ) = 0;   

};


/*! \brief Classe derivata per il modello Dirichlet-Categorical
 *
 *  Classe che rappresenta un documento per il problema del topic modeling.
 *  Gestisce le parole e si occupa del campionamento delle etichette per il topic a cui assegnare ogni parola.
 *  Gestisce i parametri del modello specifici del documento: \f$ \alpha, \pi_{j}, \bar(\pi)_{jl}, \bar(\pi)_{jr}, m_{j}, \bar{m}_{jl}, \bar{m}_{jl} \f$; 
 *  si occupa del campionamento di queste quantità.
 *  Tiene traccia dei conteggi delle parole nei topic
 * \date Febbraio 2016
 */


template <unsigned int DIM=1>
class CategoricalDocument final: public GenericDocument <TypeCategorical<DIM>,DIM>{

   	public:
	
		 /*! \brief Statistiche per aggiornare gli iperparametri della distribuzione del parametro latente
         */
		using STAT  = TypeCategorical<1>::STAT;
		
		/*! \brief Parametro latente, vettore di pesi delle parole distinte nel cluster
         */
		using THETA = TypeCategorical<1>::THETA;
		
		/*! \brief Singolo dato
        */
		using POINT = TypeCategorical<1>::Point;
		
		/*! \brief Identificativo del cluster
        */
		using ClusterID = unsigned int;
	
    private:

		/*! \brief Contenitore dei dati.
		 *  Per ogni cluster ho il subcluster sinistro e destro: nella mappa la chiave é il dato e il valore mappato é il numero di volte che il dato compare in quel documento,
		 *  in quel cluster, in quel subcuster
         */
		unordered_map < ClusterID , pair < unordered_map < POINT, unsigned int>, unordered_map < POINT, unsigned int> > >  Zeta;

	    /*! \brief Vocabolario delle parole distinte del documento
		 *  La chiave nella mappa è il dato, il valore mapparo é il numero di volte che il dato compare nel documento
		 */
		unordered_map <POINT, unsigned int> Vocabulary;
		
		/*! \brief Iperparametro del processo di Dirichlet che governa il documento
		 */
		double alpha;
		
		/*! \brief Numero di parole contenute nel documento
		 */
		unsigned int Nj;
		
		/*! \brief Vettore di pesi del cluster specifici del documento
         */		
		vector<double> Pi;
		
		/*! \brief Vettore di pesi del subcluster sinistro specifici del documento
         */	
		vector<double> PiLeft;
		
		/*! \brief Vettore di pesi del subcluster destro specifici del documento
         */	
		vector<double> PiRight;

        /*! \brief Vettore di conteggi di dimensione K: l'elemento in posizione k indica il numero di dati del documento j nel cluster k
         */		
		vector<unsigned int> WordCount;
		
		/*! \brief Vettore di conteggi di dimensione K: l'elemento in posizione k indica il numero di dati del documento j nel subcluster sinistro del cluster k
         */	
		vector<unsigned int> WordCountLeft;
		
		/*! \brief Vettore di conteggi di dimensione K: l'elemento in posizione k indica il numero di dati del documento j nel subcluster destro del cluster k
         */
		vector<unsigned int> WordCountRight;
		
		/*! \brief Vettore dei tavoli di dimensione K : l'elemento in posizione k indica il numero di tavoli nel ristorante j che servono il piatto k
		 */
		vector<unsigned int> LocalTable;
		
		/*! \brief Vettore dei tavoli di dimensione K : l'elemento in posizione k indica il numero di tavoli nel ristorante j che servono il piatto k left
		 */
		vector<unsigned int> LocalTableLeft;
		
		/*! \brief Vettore dei tavoli di dimensione K : l'elemento in posizione k indica il numero di tavoli nel ristorante j che servono il piatto k right
		 */
		vector<unsigned int> LocalTableRight;

	public:
	
		/*! \brief Costruttore alternativo
	     *  \param _alpha - iperparametro del processo di Dirichlet che governa il documento
		 */  
		CategoricalDocument(double _alpha): alpha(_alpha), Nj(0){};
		
		/*! \brief Costruttore di default
		 */
	    CategoricalDocument () = default;
		
		/*! \brief Distruttore
		 */
		~CategoricalDocument();
	    
		/*! \brief Move constructor
		 */
        CategoricalDocument(CategoricalDocument &&doc);
        
		/*! \brief Copy constructor
		 */
		CategoricalDocument (const CategoricalDocument &doc) = default;

        /*! \brief Copy assignment operator
		 */
		CategoricalDocument& operator=(const CategoricalDocument &doc);
		
        /*! \brief Move assignment operator 
		 */
		CategoricalDocument& operator=(CategoricalDocument &&doc);

        /*! \brief Aggiorna i pesi dei topic specifici del documento; si veda equazione (3.5) in Relazione_Parisi_Perego.pdf
		 *  \param _AllBeta - pesi globali dei topic
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		void UpdatePi (const vector<double>& _AllBeta, omprng& Gen );
		
		/*! \brief Aggiorna i pesi, specifici del documento, dei subtopic del topic k; ; si veda equazione (3.9) in Relazione_Parisi_Perego.pdf
		 *  \param _BetaLeft - peso globale del subtopic sinistro del topic k
		 *  \param _BetaRight - peso globale dei subtopic destro del topic k
		 *  \param k - id del topic
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		void UpdatePiSub (const double _BetaLeft, const double _BetaRight, const ClusterID k, omprng& Gen); 
		
	   	/*! \brief Aggiorna i pesi, specifici del documento, di tutti i subtopic, ; si veda equazione (3.9) in Relazione_Parisi_Perego.pdf
		 *  \param _BetaLeft - pesi globali dei subtopic sinistri
		 *  \param _BetaRight - pesi globali dei subtopic destri 
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
	    void UpdateAllPiSub(const vector<double> _BetaLeft, const vector<double> _BetaRight, omprng& Gen); 

		/*! \brief Aggiorna i tavoli; si veda equazione (3.9) in Relazione_Parisi_Perego.pdf
		 *  \param _stirling - numeri di stirling
		 *  \param _Beta - pesi globali dei topic
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
    	void UpdateLocalTable (const vector<long double>& _stirling, const vector<double>& _Beta, omprng& Gen);
		
		/*! \brief Aggiorna i tavoli dei subtopic del topic k; si veda equazione (3.12) in Relazione_Parisi_Perego.pdf
		 *  \param _stirling - numeri di stirling
		 *  \param _BetaLeft - peso globale del subtopic sinistro del topic k
		 *  \param _BetaRight - peso globale del subtopic destro del topic k
		 *  \param k - id del topic
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		void UpdateLocalTableSub_OneCluster(const vector<long double>& _stirling, const double _BetaLeft, const double _BetaRight, const ClusterID k, omprng& Gen);
       
        /*! \brief Aggiorna i tavoli dei subtopic di tutti i topic; si veda equazione (3.12) in Relazione_Parisi_Perego.pdf
		 *  \param _stirling - numeri di stirling
		 *  \param _BetaLeft - pesi globali dei subtopic sinistri
		 *  \param _BetaRight - pesi globali dei subtopic destri 
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
     	void UpdateAllLocalTableSub (const vector<long double>& _stirling, const vector<double>& _BetaLeft, const vector<double>& _BetaRight, omprng& Gen);

        /*! \brief Aggiorna l'etichetta per il topic di una parola, campionata con il metodo Sampling; si veda equazione (3.7) in Relazione_Parisi_Perego.pdf
		 *  \param _ThetaId - vettore dei pesi della parola id in tutti i topic
		 *  \param _VettId - id della parola 
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
	    void UpdateZeta(const THETA& _ThetaId, const POINT _VettId, omprng& Gen);
    
		/*! \brief Aggiorna l'etichetta per il topic e per il subtopic di una parola, campionate con il metodo Sampling
		 *  si vedano equazioni (3.7) - (3.11) in Relazione_Parisi_Perego.pdf
		 *  \param _ThetaId - vettore dei pesi della parola id in tutti i topic
		 *  \param _ThetaIdLeft - vettore dei pesi della parola id in tutti i subtopic sinistri
		 *  \param _ThetaIdRight - vettore dei pesi della parola id in tutti i subtopic destri
		 *  \param _VettId - id della parola 
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		void UpdateZeta_and_Sub(const THETA& _ThetaId, const THETA& _ThetaIdLeft,const THETA& _ThetaIdRight, const unsigned int _VettId, omprng& Gen);
	    
		/*! \brief Distribuisce la parola id nei subtopic del topic k, dopo aver campionato l'etichetta del subtopic con il metodo Sampling
		 *   si veda equazione (3.11) in Relazione_Parisi_Perego.pdf
		 *  \param _ThetaIdLeft - vettore dei pesi della parola id in tutti i subtopic sinistri
		 *  \param _ThetaIdRight - vettore dei pesi della parola id in tutti i subtopic destri
		 *  \param id - id della parola 
		 *  \param nidjk - numero di volte che la parola id nel documento j è capitata nel topic k
		 *  \param k - id del topic
		 *  \param Gen - generatore di numeri casuali in parallelo
		 */
		void UpdateZetaSub(const THETA& _ThetaIdLeft,const THETA& _ThetaIdRight ,const POINT id, const unsigned int nidjk, const ClusterID k, omprng& Gen);

		/*! \brief Metodo che serve nelle mosse di M-H per rimuovere il topic k
		 *  \param _k - id del topic
		 */
		void UpdateZeta(const ClusterID _k);
		
		/*! \brief Metodo che serve nelle mosse di M-H per rimuovere due topic
		 *  \param _k1 - id del topic
		 *  \param _k2 - id del topic
		 */
		void UpdateZeta(const ClusterID _k1, const ClusterID _k2);

		/*! \brief Visualizza il numero di parole nel documento j
		 *  \return Numero di parole nel documento j
		 */
		unsigned int ViewNj() const ;
		
        /*! Estrae gli identificativi delle parole nel documento j
		 *  \param _VettId - viene riempito con gli id delle parole nel documento j
		 */
		void ViewData( vector<POINT>& _VettId) const;
		
		/*! \brief Visualizza i conteggi necessari per aggiornare i parametri latenti dei subtopic
		 *  \param _k - id del topic
		 *  \param _counts4cleft - conteggi per parametro latente del subtopic sinistro
		 *  \param _counts4cright - conteggi per parametro latente del subtopic destro
		 */
		void ViewCounts4c(ClusterID _k, STAT& _counts4cleft, STAT& _counts4cright);


	    /*! \brief Estrae il numero di tavoli in uno specifico topic
		 *  \param _k - id del topic
		 *  \return Numero di tavoli nel ristorante j che servono il piatto k
		 */
		unsigned int ViewNumTableID(const ClusterID _k) const;
		
	    /*! \brief Estrae il numero di tavoli del subtopic sinistro di uno specifico topic
		 *  \param _k - id del topic
		 *  \return Numero di tavoli nel ristorante j che servono il piatto k left
		 */
        unsigned int ViewNumTableLeftID(const ClusterID _k) const;
	
	    /*! \brief Estrae il numero di tavoli del subtopic destro di uno specifico topic
		 *  \param _k - id del topic
		 *  \return Numero di tavoli nel ristorante j che servono il piatto k destro
		 */
		unsigned int ViewNumTableRightID(const ClusterID _k) const;

	    /*! \brief Estrae il numero di parole nel topic k
		 *  \param _k - id del topic
		 *  \return Numero di parole nel topic k
		 */
	    unsigned int ViewDataCountID(const ClusterID _k) const;
		
	     /*! \brief Estrae il numero di parole nel subtopic sinistro del topic k
		 *  \param _k - id del topic
		 *  \return Numero di parole nel subtopic sinistro del topic k
		 */
		unsigned int ViewDataCountLeftID(const ClusterID _k)const;
		
		/*! \brief Estrae il numero di parole nel subtopic destro del topic k
		 *  \param _k - id del topic
		 *  \return Numero di paorle nel subtopic destro del topic k
		 */
		unsigned int ViewDataCountRightID(const ClusterID _k)const;
		
		/*! \brief Azzera il conteggio delle parole nel topic k
		 *  \param k - id del topic
		 */
		void ResetDataCountSub(const ClusterID k);

		/*! \brief Estrae il vettore del numero di parole in ogni topic
		 *  \param _WordCount - vettore del numero di parole in ogni topic
		 */
		void ViewDataCount(vector<unsigned int>& _WordCount) const;
		
		
		/*! \brief Estrae il vettore del numero di parole in ogni subtopic sinistro
		 *  \param _WordCountLeft - vettore del numero di parole in ogni subtopic sinistro
		 */
		void ViewDataCountLeft(vector<unsigned int>& _WordCountLeft) const;
		
		
		/*! \brief Estrae il vettore del numero di parole in ogni subtopic destro
		 *  \param _WordCountRight - vettore del numero di parole in ogni subtopic destro
		 */
		void ViewDataCountRight(vector<unsigned int>& _WordCountRight) const;
		
		/*! \brief Estrae identificativi e conteggi delle parole nel topic k
		 *  \param _nidjk - struttura che contiene identificativi e conteggi delle parole nel topic k
		 *  \param _k - id del topic
		 */
		void ViewIdCounts(vector<pair<POINT,unsigned int>>& _nidjk,const ClusterID _k);

		/*! \brief Estrae il topic k
		 *  \param _k - id del topic da estrarre
		 *  \param _Cluster - struttura in cui viene estratto il topic
		 */ 
		void ViewCluster(const ClusterID _k, pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>>& _Cluster) ;        

		/*! \brief Estrae il peso, specifico del documento, del topic k
		 *  \param _k - id del topic
		 *  \return Peso specifico del documento del topic k
		 */
		double ViewPiID(const ClusterID _k) const;
		
		/*! \brief Estrae il peso, specifico del documento, del subtopic sinistro del topic k
		 *  \param _k - id del topic
		 *  \return Peso specifico del documento del subtopic sinistro del topic k
		 */
		double ViewPiLeftID(const ClusterID _k) const;
		
		/*! \brief Estrae il peso, specifico del documento, del subtopic destro del topic k
		 *  \param _k - id del topic
		 *  \return Peso specifico del documento del subtopic destro del topic k
		 */
		double ViewPiRightID(const ClusterID _k) const;		
		
		/*! \brief Estrae il vettore di pesi dei topic specifici dei documenti
		 *  \param _pi - vettore di pesi dei topic specifici dei documenti
		 */
		void ViewPi( vector<double>& _pi) const;
		
		/*! \brief Estrae il vettore di pesi dei subotpic sinistri specifici dei documenti
		 *  \param _pi - vettore di pesi dei subotopic sinistri specifici dei documenti
		 */
		void ViewPiLeft( vector<double>& _pi) const;
		
		/*! \brief Estrae il vettore di pesi dei subotpic destri specifici dei documenti
		 *  \param _pi - vettore di pesi dei subotopic destri specifici dei documenti
		 */
		void ViewPiRight( vector<double>& _pi) const;
		
		/*! \brief Imposta il parametro di concentrazione del processo di Dirichlet che governa il documento
		 *  \param _alpha - parametro di concentrazione del processo di Dirichlet che governa il documento
		 */
		void SetAlpha (const double _alpha);
		
		/*! \brief Imposta il numero di parole nel documento j
		 *  \param _Nj - numero di parole nel documento j
		 */
		void SetNj (const unsigned int _Nj);
		
		/*! \brief Imposta il vettore di pesi dei topic specifici dei documenti
		 */
		void SetPi(vector<double>& _Pi);

        /*! \brief Inserisce un nuovo topic 
		 *  \param NewCluster - il nuovo topic
		 *  \param _Pi - peso del nuovo topic specifico del documento
		 *  \param _PiLeft - peso del subtopic sinistro nuovo topic, specifico del documento
		 *  \param _PiRight - peso del subtopic destro del nuovo topic, specifico del documento
		 *  \param _WordCount - numero di parole nel nuovo topic
		 *  \param _WordCountLeft - numero di parole nel subtopic sinistro del nuovo topic
		 *  \param _WordCountRight - numero di parole nel subtopic destro del nuovo topic
		 *  \param _LocalTable - numero di tavoli che servono il nuovo piatto nel ristorante j
		 *  \param _LocalTableLeft - numero di tavoli che servono il nuovo piatto left nel ristorante j
		 *  \param _LocalTable - numero di tavoli che servono il nuovo piatto right nel ristorante j
		 */
		void InsertNewCluster(const pair < unordered_map < POINT, unsigned int>, unordered_map < POINT, unsigned int> >& NewCluster, const double _Pi, const double _PiLeft, 
	                      const double _PiRight, const unsigned int _WordCount, const unsigned int _WordCountLeft, const unsigned int _WordCountRight, const unsigned int _LocalTable,
						  const unsigned int _LocalTableLeft, const unsigned int _LocalTableRight);
						  
		/*! \brief Rimuove i topic con identificativo presente nel vettore in ingresso
		 *  \param _k - vettore con gli id dei topic da eliminare
		 */
		void RemoveCluster(const vector<ClusterID>& _k);
		
		/*! \brief Rimuove un topic
		 *  \param _k - id del topic da eliminare
		 */
		void RemoveCluster(const  ClusterID _k);
		
		/*! \brief Rimuove due topic
		 *  \param _k1 - id del topic da eliminare
		 *  \param _k2 - id del topic da eliminare
		 */
		void RemoveCluster(ClusterID _k1, ClusterID _k2);

        /*! \brief Verifica se un topic ha il subtopic sinistro vuoto
		 *  \param _k - id del topic di cui controllare il subtopic
		 */
		unsigned int CheckLeftSubcluster(const ClusterID _k);
		
		/*! \brief Verifica se un topic ha il subtopic destro vuoto
		 *  \param _k - id del topic di cui controllare il subtopic
		 */
		unsigned int CheckRightSubcluster(const ClusterID _k);

		/*! \brief Estrae le etichette associate alle parole
		 *  \param Data - struttura in cui estrarre le etichette
		 */
		void ViewLabel(vector<pair<POINT,ClusterID>>& Data);
	
		//void PrintData();

	//da togliere
	//inline void ViewDocument();
	
	    /*! \brief Acquisisce i dati
		 *  \param SSTR - contiene id della parola e numero di volte che parola compare nel documento
		 */
		void SetDataset(std::istringstream& SSTR);
		
		/*! \brief Smista le parole nel contenitore Zeta
		 *  \param _K - numero iniziale di topic
		 *  \param Gen - generatore di numeri casuali
		 */
		unsigned int SortData(unsigned int _K, omprng& Gen);


    private:
		
		/*! \brief Aggiorna i conteggi delle parole nei topic
		 */
		void UpdateDataCount();
		
		/*! \brief Campionamento da distribuzione categorica per l'etichetta del topic o del subtopic
		 *  \param _temp_counts - vettore che contiene i conteggi della parola id nei topic
		 *  \param _Weights - pesi con cui campionare le etichette
		 *  \param _nidj - numero di volte che la parola id compare nel documento j
		 *  \param Gen - generatore di numeri casuali
		 */
		void Sampling (std::vector<unsigned int>& _temp_counts, std::vector<double>& _Weights, unsigned int _nidj, omprng& Gen);

};

//Acquisisce i dati
template<unsigned int DIM> void CategoricalDocument<DIM>::SetDataset(std::istringstream& SSTR){

	POINT WordId;
	unsigned int nidj;
	
	SSTR >> WordId;
	SSTR >> nidj;
	
	Vocabulary.insert(std::make_pair(WordId,nidj));
	
	Nj += nidj;	
}

//Smista i dati nel contenitore Zeta
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::SortData(unsigned int _K, omprng& Gen){
	
	Gen.setNumThreads(1);

    pair<unordered_map <POINT, unsigned int>,unordered_map <POINT, unsigned int>> Cluster_k;

//inserisco in Zeta K cluster vuoti
    for(size_t k=0; k<_K; ++k)
       Zeta.insert({k,Cluster_k});

//devo scorrere il Vocabolario. Per ogni id, devo estrarre un etichetta da 0 a K-1 in modo uniforme, per tutte le volte che l'id compare, cioè per nid volte

    unsigned int label;
    unsigned int sub_label;

    for (auto it = Vocabulary.cbegin(); it != Vocabulary.cend(); ++it){

	    for(size_t n=0; n< it->second; ++n){

	       label = Gen.runifdiscrete(_K);
	       sub_label = Gen.rbernoulli(1.0/2);

	        if(sub_label == 1){  //inserisco a sinistra, a destra metto 0
	         if(!Zeta[label].first.insert(std::make_pair(it->first,1)).second) //se inserimento non ha successo, id già presente
		        Zeta[label].first[it->first] += 1;
			 else // se inserisco nuovo id, va anche a dx
			    Zeta[label].second.insert(std::make_pair(it->first,0));
	        }
	        else{		// viceversa
	        if(!Zeta[label].second.insert(std::make_pair(it->first,1)).second) //se inserimento non ha successo, id già presente
		    Zeta[label].second[it->first] += 1;
			else // se inserisco nuovo id, va anche a sx
			 Zeta[label].first.insert(std::make_pair(it->first,0));
		    }
	    }
    }

    LocalTable.reserve(_K);
    LocalTableLeft.reserve(_K);
    LocalTableRight.reserve(_K);

    for(size_t i=0; i<_K; ++i){
       LocalTable.push_back(0);
       LocalTableLeft.push_back(0);
       LocalTableRight.push_back(0);
    }

    this-> UpdateDataCount();

	return Nj;
}

// Distruttore
template<unsigned int DIM> CategoricalDocument<DIM>::~CategoricalDocument(){}   


//Aggiorna i pesi dei topic specifici del documento
template<unsigned int DIM> void CategoricalDocument<DIM>::UpdatePi(const vector<double>& _AllBeta, omprng& Gen ){

    vector<double> dir_sampled(_AllBeta.size());
    vector<double> params;
	params.reserve(_AllBeta.size());

	for( size_t i=0; i< (_AllBeta.size()-1); ++i)
	    params.push_back(alpha*_AllBeta[i]+ WordCount[i]);
    params.push_back(alpha*_AllBeta[_AllBeta.size()-1]);

	Gen.rdirichlet(params,dir_sampled);

	Pi = dir_sampled;

}

 // Move constructor
template<unsigned int DIM> CategoricalDocument<DIM>::CategoricalDocument(CategoricalDocument<DIM> &&doc):
	Nj(doc.Nj)
	,alpha(doc.alpha)
	,Zeta(doc.Zeta)
	,Vocabulary(doc.Vocabulary)
	,Pi(doc.Pi)
	,PiLeft(doc.PiLeft)
	,PiRight(doc.PiRight)
	,WordCount(doc.WordCount)
	,WordCountLeft(doc.WordCountLeft)
	,WordCountRight(doc.WordCountRight)
	,LocalTable(doc.LocalTable)
	,LocalTableLeft(doc.LocalTableLeft)
	,LocalTableRight(doc.LocalTableRight)
	{}


// Operatore di assegnamento
template<unsigned int DIM> CategoricalDocument<DIM>& CategoricalDocument<DIM>::operator=(const CategoricalDocument<DIM>& doc){

    Zeta = doc.Zeta;
    Vocabulary = doc.Vocabulary;
    Nj = doc.Nj;
    Pi = doc.Pi;
    PiLeft = doc.PiLeft;
	PiRight = doc.PiRight;
	WordCount = doc.WordCount;
	WordCountLeft = doc.WordCountLeft;
	WordCountRight = doc.WordCountRight;
	LocalTable = doc.LocalTable;
	LocalTableLeft = doc.LocalTableLeft;
	LocalTableRight = doc.LocalTableRight;
	alpha = doc.alpha;

	return *this;
}

// Move assignment operator
template<unsigned int DIM> CategoricalDocument<DIM>& CategoricalDocument<DIM>::operator=(CategoricalDocument<DIM>&& doc) {

    Zeta = doc.Zeta;
    Vocabulary = doc.Vocabulary;
    Nj = doc.Nj;
	Pi = doc.Pi;
	PiLeft = doc.PiLeft;
	PiRight = doc.PiRight;
	WordCount = doc.WordCount;
	WordCountLeft = doc.WordCountLeft;
	WordCountRight = doc.WordCountRight;
	LocalTable = doc.LocalTable;
	LocalTableLeft = doc.LocalTableLeft;
	LocalTableRight = doc.LocalTableRight;
	alpha = doc.alpha;

	return *this;
}

//Aggiorna i pesi, specifici del documento, dei subtopic del topic k
template<unsigned int DIM> void CategoricalDocument<DIM>::UpdatePiSub (const double _BetaLeft, const double _BetaRight, const ClusterID k, omprng& Gen){

	double first_param= alpha * _BetaLeft+ WordCountLeft[k];
	double second_param= alpha * _BetaRight + WordCountRight[k];

	vector<double> params{first_param,second_param};
	vector<double> dirichlet_sampled(2);

    Gen.rdirichlet(params,dirichlet_sampled);

	PiLeft[k] = dirichlet_sampled[0];
	PiRight[k] = dirichlet_sampled[1];

}

//Aggiorna i pesi, specifici del documento, di tutti i subtopic
template<unsigned int DIM> void CategoricalDocument<DIM>::UpdateAllPiSub(const vector<double> _BetaLeft, const vector<double> _BetaRight, omprng& Gen){

    if(_BetaLeft.size() != _BetaRight.size()) {
	   std::cerr<< " In UpdatePiSub, BetaLeft e BetaRight must have the same dimension"<<std::endl;
	   exit(1);
	}

    unsigned int K= _BetaLeft.size();

	PiLeft.clear();
	PiLeft.resize(K);
	PiRight.clear();
	PiRight.resize(K);

	for(size_t k=0; k<K; ++k)
        this -> UpdatePiSub(_BetaLeft[k],_BetaRight[k],k,Gen);

}


//Aggiorna i conteggi delle parole nei topic
template<unsigned int DIM> void CategoricalDocument<DIM>::UpdateDataCount(){

    vector<unsigned int> WordCount_new;
    vector<unsigned int> WordCountLeft_new;
	vector<unsigned int> WordCountRight_new;

	WordCount_new.resize(Zeta.size());
	WordCountLeft_new.resize(Zeta.size());
	WordCountRight_new.resize(Zeta.size());

	for(auto iter_out= Zeta.cbegin(); iter_out != Zeta.cend(); ++iter_out){

        auto first_map= iter_out->second.first;
	    auto second_map= iter_out->second.second;
	    unsigned int first_sum=0;
		unsigned int second_sum=0;

		for(auto iter_in= first_map.cbegin(); iter_in != first_map.cend(); ++iter_in)
			first_sum += iter_in->second;
		WordCountLeft_new[iter_out->first]=first_sum;

		for(auto iter_in= second_map.cbegin(); iter_in != second_map.cend(); ++iter_in)
			second_sum += iter_in->second;
		WordCountRight_new[iter_out->first]=second_sum;

		WordCount_new[iter_out->first]=first_sum + second_sum; //njk

	}

	WordCount=WordCount_new;
    WordCountLeft=WordCountLeft_new;
	WordCountRight=WordCountRight_new;

}

//Azzera il conteggio delle parole nel topic k
template<unsigned int DIM> void CategoricalDocument<DIM>::ResetDataCountSub(const ClusterID k){
    WordCountLeft[k] = 0;
	WordCountRight[k] = 0;
}

//Aggiorna i tavoli 
template<unsigned int DIM> void CategoricalDocument<DIM>::UpdateLocalTable( const vector<long double>& _stirling, const vector<double>& _Beta, omprng& Gen){

    unsigned int K= _Beta.size();
	unsigned int mjk=0;
	LocalTable.clear();
	LocalTable.reserve(K);

	for(size_t i=0; i<K; ++i){
	    mjk= Antoniak(alpha,_Beta[i],WordCount[i],_stirling,Gen);
		LocalTable.push_back(mjk);
	}

}

//Aggiorna i tavoli dei subtopic del topic k
template<unsigned int DIM>
void CategoricalDocument<DIM>::UpdateLocalTableSub_OneCluster(const vector<long double>& _stirling, const double _BetaLeft, const double _BetaRight, const ClusterID k, omprng& Gen){

    if(k > WordCountLeft.size()) {
	  std::cerr<< " In UpdateLocalTableSub, you have selected a non-existing cluster "<<std::endl;
	  exit(1);
	}

    unsigned int mjkl=0;
    unsigned int mjkr=0;

    mjkl=Antoniak(alpha,_BetaLeft,WordCountLeft[k],_stirling, Gen);
    LocalTableLeft[k] = mjkl;

    mjkr=Antoniak(alpha,_BetaRight,WordCountRight[k],_stirling, Gen);
    LocalTableRight[k] = mjkr;

}

// Aggiorna i tavoli dei subtopic di tutti i topic
template<unsigned int DIM>
void CategoricalDocument<DIM>::UpdateAllLocalTableSub(const vector<long double>& _stirling, const vector<double>& _BetaLeft, const vector<double>& _BetaRight, omprng& Gen){

    if(_BetaLeft.size() != _BetaRight.size()) {
	   std::cerr<< " In UpdateLocalTableSub, BetaLeft e BetaRight must have the same dimension"<<std::endl;
	   exit(1);
	}

	unsigned int K = _BetaLeft.size();

	LocalTableLeft.clear();
	LocalTableLeft.resize(K);
	LocalTableRight.clear();
	LocalTableRight.resize(K);

	for(size_t i=0; i<K;++i)
       this->UpdateLocalTableSub_OneCluster(_stirling, _BetaLeft[i],_BetaRight[i],i,Gen);

}

//Visualizza il numero di parole nel documento j
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::ViewNj() const {
	return Nj;
}


// Estrae gli identificativi delle parole nel documento j
template<unsigned int DIM> void CategoricalDocument<DIM>::ViewData(vector<POINT>& _VettId) const{

    if(!(_VettId.empty())){
	    std::cerr<<"You must pass an empty vector (memory NOT reserved) to  Document::ViewData"<<std::endl;
		exit(1);
	}

    unsigned int v= Vocabulary.size();
	_VettId.reserve(v);

    for(auto it= Vocabulary.cbegin(); it!= Vocabulary.cend(); ++it)
	    _VettId.push_back(it->first);
}


//Aggiorna l'etichetta per il topic di una parola
template<unsigned int DIM> void CategoricalDocument<DIM>::UpdateZeta(const THETA& _ThetaId, const unsigned int _VettId, omprng& Gen){

	vector<unsigned int> temp_counts;
	unsigned int K= _ThetaId.size();
	vector<double> Weights;
	Weights.reserve(K);

	if(Pi.empty()){
	  std::cerr<< "Pi has not been allocated , if you proceed Segmentation Fault"<<std::endl;
	  exit(1);
	}

	if(PiLeft.empty()){
	  std::cerr<< "PiLeft has not been allocated , if you proceed Segmentation Fault"<<std::endl;
	  exit(1);
	}

	for( vector<double>::size_type i=0; i<K;++i){
		Weights.push_back( Pi[i]*_ThetaId[i]);
	}

	double weight_sum=0.0;
    weight_sum= Kahan_algorithm(Weights);

    for(auto& i: Weights)
        i= i/weight_sum;

    double weight_sum_after=0.0;
    weight_sum_after = Kahan_algorithm(Weights);

    if( !((std::fabs(1.0 - weight_sum_after)/std::fabs(1.0)) < 9.99e-10)){
       std::cerr<<"Error in UpdateZeta"<<std::endl;
       std::cerr<<"weight_sum_after is not equal to 1 "<<std::endl;
       exit(1);
    }

	unsigned int nidj = Vocabulary[_VettId];

	this -> Sampling (temp_counts,Weights,nidj,Gen);

	if(temp_counts.empty()){
	  std::cerr<<"Error in UpdateZeta, Categorical function has not sampled temp_counts"<<std::endl;
	  exit(1);
	}

	for( size_t i=0; i<K;i++){
        Zeta[i].first[_VettId]= temp_counts[i];
		Zeta[i].second[_VettId]= 0;

	}

	this-> UpdateDataCount();
}

//Aggiorna l'etichetta per il topic e per il subtopic di una parola
template<unsigned int DIM>
void CategoricalDocument<DIM>::UpdateZeta_and_Sub(const THETA& _ThetaId, const THETA& _ThetaIdLeft, const THETA& _ThetaIdRight, const unsigned int _VettId, omprng& Gen){

	if((_ThetaId.size()!= _ThetaIdLeft.size()) || (_ThetaId.size()!=_ThetaIdRight.size())){
	    std::cerr<<"ThetaId, Left e Right in UpdateZeta must have the same dimension"<<std::endl;
		exit(1);
	}

	vector<unsigned int> temp_counts;
	vector<unsigned int> temp_counts_left;
	unsigned int temp_counts_right=0;
	unsigned int K= _ThetaId.size();

	vector<double> WeightsZeta;
	WeightsZeta.reserve(K);

	vector<double> WeightsZetaSub;
	WeightsZetaSub.resize(1);

	double pi_theta_left = 0.0;
	double pi_theta_right = 0.0;

	if(Pi.empty()){
		std::cerr<< "Pi has not been allocated , if you proceed Segmentation Fault"<<std::endl;
		exit(1);
	}

	if(PiLeft.empty()){
	    std::cerr<< "PiLeft has not been allocated , if you proceed Segmentation Fault"<<std::endl;
		exit(1);
	}


	for( vector<double>::size_type i=0; i<K;++i){
		WeightsZeta.push_back( Pi[i]*_ThetaId[i]);
	}

	double weight_sum=0.0;
    weight_sum= Kahan_algorithm(WeightsZeta);

    for(auto& i: WeightsZeta)
        i= i/weight_sum;

    double weight_sum_after=0.0;
    weight_sum_after = Kahan_algorithm(WeightsZeta);

    if( !((std::fabs(1.0 - weight_sum_after)/std::fabs(1.0)) < 9.99e-10)){
        std::cerr<<"Errore in Update Zeta and SUb"<<std::endl;
        exit(1);
    }

	unsigned int nidj = Vocabulary[_VettId];

	this -> Sampling (temp_counts,WeightsZeta,nidj,Gen);

	if(temp_counts.empty()){
	    std::cerr<<"Error in UpdateZeta_andSub, Categorical function hasn't sampled temp_counts"<<std::endl;
		exit(1);
	}

	for(size_t i=0; i<K;i++){

        pi_theta_left = PiLeft[i]*_ThetaIdLeft[i];
		pi_theta_right = PiRight[i]*_ThetaIdRight[i];

		weight_sum = pi_theta_left + pi_theta_right;

		WeightsZetaSub[0]= pi_theta_left / weight_sum;

		this -> Sampling(temp_counts_left,WeightsZetaSub,temp_counts[i],Gen);

		if(temp_counts_left.empty()){
		   std::cerr<<"Error in UpdateZeta_andSub, Categorical function hasn't sampled temp_counts_left" <<std::endl;
		   exit(1);
		}

		temp_counts_right = temp_counts[i]- temp_counts_left[0];

		Zeta[i].first[_VettId]= temp_counts_left[0];
		Zeta[i].second[_VettId]= temp_counts_right;
		temp_counts_left.clear();

	}

	this-> UpdateDataCount();
}

// questo metodo viene usato in Gibbs Subtopic (almeno per le mosse locali)
template<unsigned int DIM>
void CategoricalDocument<DIM>::UpdateZetaSub(const THETA& _ThetaIdLeft,const THETA& _ThetaIdRight ,const POINT id, const unsigned int nidjk, const ClusterID k, omprng& Gen){

	if(_ThetaIdLeft.size()!= 1){
	   std::cerr<<"In UpdateZetaSub ThetaLeft must be 1-D"<<std::endl;
	   exit(1);
	}

    vector<unsigned int> temp_counts_sub;

	vector<double> weights{PiLeft[k]*_ThetaIdLeft[0],PiRight[k]*_ThetaIdRight[0]};

	double weight_sum=0.0;

	weight_sum= Kahan_algorithm(weights);

    for (auto& i: weights)
        i= i/weight_sum;

	this -> Sampling(temp_counts_sub,weights,nidjk,Gen);

	if(temp_counts_sub.empty()){
	  std::cerr<<"Error in UpdateZetaSub, Categorical function hasn't sampled temp_counts_sub" <<std::endl;
	  exit(1);
	}

    if( ((temp_counts_sub[0] + temp_counts_sub[1] )!= nidjk )){
	  std::cerr<< "Error in UpdateZetaSub, counts do not match nidjk"<<std::endl;
	  exit(1);
	}

	Zeta[k].first[id]= temp_counts_sub[0];
	Zeta[k].second[id]= temp_counts_sub[1];

	// riguardo word count devo aggiornare solo WCL e WCR del cluster proposto che sto considerando
	WordCountLeft[k] += temp_counts_sub[0];
	WordCountRight[k] += temp_counts_sub[1];
 }


//Visualizza i conteggi necessari per aggiornare i parametri latenti dei subtopic
template<unsigned int DIM>
void CategoricalDocument<DIM>::ViewCounts4c( const ClusterID _k, STAT& _counts4cleft, STAT& _counts4cright ) {
    
	if(_counts4cleft.size() != _counts4cright.size() ){
		std::cerr<<" counts4cleft e counts4cright in ViewCounts4c must have the same dimension"<<std::endl;
		exit(1);
	}
	
	unsigned int W = _counts4cleft.size();

	_counts4cleft.assign(W,0);
	_counts4cright.assign(W,0);

	auto left_map= Zeta[_k].first;
	auto right_map= Zeta[_k].second;

	for(auto it= left_map.cbegin(); it!=left_map.cend(); ++it )
		_counts4cleft[(it->first)-1]=it->second;

	for(auto it= right_map.cbegin(); it!=right_map.cend(); ++it )
		_counts4cright[(it->first)-1]=it->second;
}

//Metodo che serve nelle mosse di M-H per rimuovere il topic k
template<unsigned int DIM> void CategoricalDocument<DIM>::UpdateZeta(const ClusterID _k){

	unordered_map < ClusterID , pair < unordered_map < POINT, unsigned int>, unordered_map < POINT, unsigned int> > >  ZetaNew;

	//inserisco gli elementi prima di quello da togliere (stesse chiavi)
	for(size_t i=0; i < _k ; ++i)
	    ZetaNew.insert(std::make_pair(i,Zeta[i]));
	//inserisco gli elementi dopo quello da togliere (chiave scalate all'indietro di 1)
	for(size_t i= (_k+1); i < Zeta.size() ; ++i)
	    ZetaNew.insert(std::make_pair(i-1,Zeta[i]));

    Zeta.swap(ZetaNew); //non c'è più bisogno di fare erase, dopo lo swap la vecchia zeta viene distrutta
}

//Metodo che serve nelle mosse di M-H per rimuovere due topic
template<unsigned int DIM> void CategoricalDocument<DIM>::UpdateZeta(const ClusterID _k1, const ClusterID _k2){

	unordered_map < ClusterID , pair < unordered_map < POINT, unsigned int>, unordered_map < POINT, unsigned int> > >  ZetaNew;
    //inserisco gli elementi prima del primo da togliere (stesse chiavi)
	for(size_t i=0; i < _k1 ; ++i)
	    ZetaNew.insert(std::make_pair(i,Zeta[i]));
	//inserisco gli elementi tra i due da togliere (chiavi scalate in modo opportuno)
	for(size_t i=_k1+1; i < _k2 ; ++i)
	    ZetaNew.insert(std::make_pair(i-1,Zeta[i]));
	//inserisco gli elementi dopo del secondo da togliere (chiave scalate all'indietro di 1)
	for(size_t i= (_k2+1); i < Zeta.size() ; ++i)
	    ZetaNew.insert(std::make_pair(i-2,Zeta[i]));

	Zeta.swap(ZetaNew);
}

//Estrae il numero di tavoli in uno specifico topic
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::ViewNumTableID(const ClusterID _k) const{
    return LocalTable[_k];
}

//Estrae il numero di tavoli del subtopic sinistro di uno specifico topic
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::ViewNumTableLeftID(const ClusterID _k) const{
    return LocalTableLeft[_k];
}

//Estrae il numero di tavoli del subtopic destro di uno specifico topic
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::ViewNumTableRightID(const ClusterID _k) const{
    return LocalTableRight[_k];
}

//Estrae il vettore del numero di parole in ogni topic
template<unsigned int DIM> void CategoricalDocument<DIM>::ViewDataCount(vector<unsigned int>& _WordCount)const{
	_WordCount = WordCount;
}

//Estrae il vettore del numero di parole in ogni subtopic sinistro
template<unsigned int DIM> void CategoricalDocument<DIM>::ViewDataCountLeft(vector<unsigned int>& _WordCountLeft)const{
	_WordCountLeft = WordCountLeft;
}

//Estrae il vettore del numero di parole in ogni subtopic destro
template<unsigned int DIM> void CategoricalDocument<DIM>::ViewDataCountRight(vector<unsigned int>& _WordCountRight) const{
	_WordCountRight = WordCountRight;
}

//Estrae il numero di parole nel topic k
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::ViewDataCountID (const ClusterID _k) const{
    return WordCount[_k];
}

//Estrae il numero di parole nel subtopic sinistro del topic k
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::ViewDataCountLeftID(const ClusterID _k) const{
    return WordCountLeft[_k];
}

//Estrae il numero di parole nel subtopic destro del topic k
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::ViewDataCountRightID(const ClusterID _k)const{
    return WordCountRight[_k];
}

//Estrae il peso, specifico del documento, del topic k
template<unsigned int DIM> double CategoricalDocument<DIM>::ViewPiID(const ClusterID _k)const{
    return Pi[_k];
}

//Estrae il peso, specifico del documento, del subtopic sinistro del topic k
template<unsigned int DIM> double CategoricalDocument<DIM>::ViewPiLeftID(const ClusterID _k)const{
    return PiLeft[_k];
}

//Estrae il peso, specifico del documento, del subtopic destro del topic k
template<unsigned int DIM> double CategoricalDocument<DIM>::ViewPiRightID(const ClusterID _k)const{
    return PiRight[_k];
}

//Estrae il vettore di pesi dei topic specifici dei documenti
template<unsigned int DIM> void CategoricalDocument<DIM>::ViewPi(vector<double>& _pik) const{
    _pik = Pi;
}

//Estrae il vettore di pesi dei subotpic sinistri specifici dei documenti
template<unsigned int DIM> void CategoricalDocument<DIM>::ViewPiLeft(vector<double>& _pik) const{
    _pik = PiLeft;
}

//Estrae il vettore di pesi dei subotpic destri specifici dei documenti
template<unsigned int DIM> void CategoricalDocument<DIM>::ViewPiRight(vector<double>& _pik) const{
    _pik = PiRight;
}

//Estrae il topic k
template<unsigned int DIM>
void CategoricalDocument<DIM>::ViewCluster(const ClusterID _k, pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>>& _Cluster) {
	 _Cluster = Zeta[_k];
}

//Imposta il parametro di concentrazione del processo di Dirichlet che governa il documento
template<unsigned int DIM> void CategoricalDocument<DIM>::SetAlpha(double _alpha){
    alpha = _alpha;
}

//Imposta il numero di parole nel documento j
template<unsigned int DIM> void CategoricalDocument<DIM>::SetNj(unsigned int _Nj){
    Nj = _Nj;
}

//Imposta il vettore di pesi dei topic specifici dei documenti
template<unsigned int DIM> void CategoricalDocument<DIM>::SetPi(vector<double>& _Pi){
    Pi = _Pi;
}

//Rimuove i topic con identificativo presente nel vettore in ingresso
template<unsigned int DIM> void CategoricalDocument<DIM>::RemoveCluster(const vector<ClusterID>& _k){


    if(_k.empty())
       return;

    unordered_map < ClusterID , pair < unordered_map < POINT, unsigned int>, unordered_map < POINT, unsigned int> > >  ZetaNew;

    vector<double> PiNew;
	vector<double> PiLeftNew;
	vector<double> PiRightNew;

	double TempPi_EmptyCluster;

	if(!Pi.empty())
        TempPi_EmptyCluster = *(Pi.rbegin());

		// njk, nr parole che nel documento j è assegnato al cluster k, mi servono metodi che fanno le somme in Zeta sugli n_id_j_k_h
		//è un vettore lungo k+1, dove nj_k+1 è sempre=0
	vector<unsigned int> WordCountNew;
	vector<unsigned int> WordCountLeftNew;
	vector<unsigned int> WordCountRightNew;
		//mjk per ogni k, vengono campionati salla Antoniak
	vector<unsigned int> LocalTableNew;
	vector<unsigned int> LocalTableLeftNew;
	vector<unsigned int> LocalTableRightNew;

	bool check = true;
	ClusterID NewKey = 0;

    for(size_t i = 0; i < Zeta.size(); i++){
        for(auto it = _k.cbegin(); it != _k.cend(); it++){

            if(*it == i){
                check = false;
                break;
            }
            else {
                check = true;
            }
        }
        if(check){
            ZetaNew.insert(std::make_pair(NewKey,Zeta[i]));
	    if(!Pi.empty()){
	    	PiNew.push_back(Pi[i]);
            PiLeftNew.push_back(PiLeft[i]);
            PiRightNew.push_back(PiRight[i]);
	    }

        WordCountNew.push_back(WordCount[i]);
        WordCountLeftNew.push_back(WordCountLeft[i]);
        WordCountRightNew.push_back(WordCountRight[i]);

        LocalTableNew.push_back(LocalTable[i]);
        LocalTableLeftNew.push_back(LocalTableLeft[i]);
        LocalTableRightNew.push_back(LocalTableRight[i]);
        NewKey++;
        }

    }

    if(!Pi.empty()){
    	PiNew.push_back(TempPi_EmptyCluster);
	    Pi.swap(PiNew);
    	PiLeft.swap(PiLeftNew);
    	PiRight.swap(PiRightNew);

    }

    Zeta.swap(ZetaNew);
    WordCount.swap(WordCountNew);
    WordCountLeft.swap(WordCountLeftNew);
    WordCountRight.swap(WordCountRightNew);
    LocalTable.swap(LocalTableNew);
    LocalTableLeft.swap(LocalTableLeftNew);
    LocalTableRight.swap(LocalTableRightNew);

}

//Rimuove un topic
template<unsigned int DIM> void CategoricalDocument<DIM>::RemoveCluster(const ClusterID _k){

	 this-> UpdateZeta(_k); // creo una nuova Zeta senza il cluster _k e scambio il contenuto con quella vecchia
	 Pi.erase(Pi.begin() + _k );  //sistemo i vettori, erase distrugge l'elemento in posizione _k e il vettore si riduce di dimensione
	 PiLeft.erase(PiLeft.begin() + _k );
	 PiRight.erase(PiRight.begin() + _k );
	 WordCount.erase(WordCount.begin() + _k );
	 WordCountLeft.erase(WordCountLeft.begin() + _k );
	 WordCountRight.erase(WordCountRight.begin() + _k );
	 LocalTable.erase(LocalTable.begin() + _k );
	 LocalTableLeft.erase(LocalTableLeft.begin() + _k );
	 LocalTableRight.erase(LocalTableRight.begin() + _k );
}

//Rimuove due topic
template<unsigned int DIM> void CategoricalDocument<DIM>::RemoveCluster( ClusterID _k1, ClusterID _k2){

	if(_k1==_k2){
	  std::cerr<< "In Document::RemoveCluster(_k1,_k2), k1 and k2 cannot be the same number"<<std::endl;
	  exit(1);

	}

	if(!(_k1<_k2)){
	   unsigned int temp = _k1;
	   _k1=_k2;
	   _k2=temp;
	}

	 this-> UpdateZeta(_k1,_k2); //sistemo le chiavi
	 Pi.erase( Pi.begin() + _k1);  //sistemo i vettori, erase distrugge l'elemento in posizione _k e il vettore si riduce di dimensione
	 Pi.erase( Pi.begin() + _k2-1);  // -1 perchè dopo il primo erase sono scalate di 1 le posizioni
	 PiLeft.erase(PiLeft.begin()+ _k1);
	 PiLeft.erase(PiLeft.begin() + _k2-1);
	 PiRight.erase(PiRight.begin()+ _k1);
	 PiRight.erase(PiRight.begin() + _k2-1);
	 WordCount.erase( WordCount.begin() + _k1);
	 WordCount.erase( WordCount.begin() + _k2-1);
	 WordCountLeft.erase( WordCountLeft.begin() + _k1);
	 WordCountLeft.erase( WordCountLeft.begin() + _k2-1);
     WordCountRight.erase( WordCountRight.begin() + _k1);
	 WordCountRight.erase( WordCountRight.begin() + _k2-1);
	 LocalTable.erase( LocalTable.begin() + _k1);
	 LocalTable.erase( LocalTable.begin() + _k2-1);
	 LocalTableLeft.erase( LocalTableLeft.begin() + _k1);
	 LocalTableLeft.erase( LocalTableLeft.begin() + _k2-1);
	 LocalTableRight.erase( LocalTableRight.begin() + _k1);
	 LocalTableRight.erase( LocalTableRight.begin() + _k2-1);

}

//Inserisce un nuovo topic 
template<unsigned int DIM>
void CategoricalDocument<DIM>::InsertNewCluster(const pair < unordered_map < POINT, unsigned int>, unordered_map < POINT, unsigned int> >& NewCluster, 
												const double _Pi, const double _PiLeft, const double _PiRight, const unsigned int _WordCount, const unsigned int _WordCountLeft, 
												const unsigned int _WordCountRight, const unsigned int _LocalTable, const unsigned int _LocalTableLeft, const unsigned int _LocalTableRight){

    ClusterID NewKey=0;
    double TempPi_EmptyCluster = 0.0;
    unsigned int Pos_EmptyCluster = 0;

    TempPi_EmptyCluster = *(Pi.rbegin());
    Pos_EmptyCluster = Pi.size() - 1;
    Pi.erase(Pi.begin() + Pos_EmptyCluster);

    if(Zeta.size() == 0){
	    if(!Zeta.insert({NewKey,NewCluster}).second){
            std::cerr<<"Errore in Document: InsertNewCluster"<<std::endl;
            exit(1);
        }
    }
    else{
    for(auto it = Zeta.begin(); it != Zeta.end(); it++){
	    if(it->first > NewKey){
		NewKey = it->first;
	    }

    }
    NewKey++;
	if(!Zeta.insert({NewKey,NewCluster}).second){
        std::cerr<<"Errore in Document: InsertNewCluster"<<std::endl;
        exit(1);
    }
    }
    Pi.push_back(_Pi);
    //reinserisco il cluster vuoto
    Pi.push_back(TempPi_EmptyCluster);
    PiLeft.push_back(_PiLeft);
    PiRight.push_back(_PiRight);
		// njk, nr parole che nel documento j è assegnato al cluster k, mi servono metodi che fanno le somme in Zeta sugli n_id_j_k_h
		//è un vettore lungo k+1, dove nj_k+1 è sempre=0
    WordCount.push_back(_WordCount);
    WordCountLeft.push_back(_WordCountLeft);
    WordCountRight.push_back(_WordCountRight);

    LocalTable.push_back(_LocalTable);
    LocalTableLeft.push_back(_LocalTableLeft);
    LocalTableRight.push_back(_LocalTableRight);

}

//Verifica se un topic ha il subtopic sinistro vuoto
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::CheckLeftSubcluster(const ClusterID _k) {

    unsigned int sum = 0;
    auto left_map = Zeta[_k].first;

	for(auto it= left_map.cbegin(); it != left_map.cend(); ++it)
		sum += it ->second;

    return sum;
}

//Verifica se un topic ha il subtopic destro vuoto
template<unsigned int DIM> unsigned int CategoricalDocument<DIM>::CheckRightSubcluster(const ClusterID _k) {

    unsigned int sum = 0;
    auto right_map = Zeta[_k].second;

	for(auto it= right_map.cbegin(); it != right_map.cend(); ++it)
		sum += it ->second;

    return sum;
}

//Estrae identificativi e conteggi delle parole nel topic k
template<unsigned int DIM> void CategoricalDocument<DIM>::ViewIdCounts(vector<pair<POINT,unsigned int>>& _nidjk,const ClusterID _k){

	 vector<pair<unsigned int,unsigned int>> _nidjk_left;
	 vector<pair<unsigned int,unsigned int>> _nidjk_right;
     auto k_left_map = Zeta[_k].first;
	 auto k_right_map = Zeta[_k].second;

	 for(auto it = k_left_map.cbegin(); it!= k_left_map.cend(); ++it)
	     _nidjk_left.push_back(std::make_pair(it->first,it->second));
	 for(auto it = k_right_map.cbegin(); it!= k_right_map.cend(); ++it)
	     _nidjk_right.push_back(std::make_pair(it->first,it->second));

	if((_nidjk_left.size()) != (_nidjk_right.size())){
	  std::cerr << " You are in HDP_MCMC::logq. Error in Document::ViewIdCounts: left and right map must have the same dimension"<<std::endl;
	  exit(1);
	}

	_nidjk.resize(_nidjk_left.size(),std::make_pair(0,0));

	auto it_left = _nidjk_left.cbegin();
	auto it_right = _nidjk_right.cbegin();

    for(vector<unsigned int>::size_type i=0; i<_nidjk_left.size(); ++i){
	    if((it_left->first)!= (it_right->first)){
		  std::cerr<< " In Document.ViewIDCounts: sx and dx maps do not contain ids in the same order"<<std::endl;
		  exit(1);
		}
	    _nidjk[i].first = it_left->first  ;
		_nidjk[i].second = it_left->second + it_right->second ;
		++it_left;
		++it_right;
	}

}

//Campionamento da distribuzione categorica per l'etichetta del topic o del subtopic
template<unsigned int DIM> void CategoricalDocument<DIM>::Sampling (std::vector<unsigned int>& _temp_counts, std::vector<double>& _Weights, unsigned int _nidj, omprng& Gen){

	double r=1;
    unsigned int v=0;

	unsigned int K= _Weights.size();
	_temp_counts.reserve(K);

    for (size_t i=0; i<K; i++){

		v= Gen.rbinomial(_nidj,_Weights[i]/r);
        _temp_counts.push_back(v);
        _nidj=_nidj-v;
        r=r-_Weights[i];
    }

}

//Estrae le etichette associate alle parole
template<unsigned int DIM>
void CategoricalDocument<DIM>::ViewLabel(vector<pair<POINT,ClusterID>>& Data) {

   unsigned int K = Pi.size() - 1;
   unsigned int nidkl;
   unsigned int nidkr;

    for(auto ItVoc = Vocabulary.begin(); ItVoc != Vocabulary.end(); ItVoc ++){
        for(size_t k=0; k<K; k++){
            nidkl = (Zeta[k].first)[ItVoc->first];
            nidkr = (Zeta[k].second)[ItVoc->first];

            for(size_t nid = 0; nid< nidkl+ nidkr; nid++){
                Data.push_back(std::make_pair(ItVoc->first,k));
            }
        }
    }

}


/*
template<unsigned int DIM>
void CategoricalDocument<DIM>::PrintData(){

    std::ofstream file ("./cpp_results/Data.bin",std::ios::binary | std::ios::app);
	if (file.fail()){
		std::cerr <<"ERROR GENERATED (DATA)!" <<std::endl;
		std::cerr <<"Cannot open Data.bin" <<std::endl;
		exit(1);
	}

    vector<POINT> Data;
    unsigned int K = Pi.size() - 1;
    unsigned int nidkl;
    unsigned int nidkr;

    for(auto ItVoc = Vocabulary.begin(); ItVoc != Vocabulary.end(); ItVoc ++){
        for(size_t k=0; k<K; k++){
            nidkl = (Zeta[k].first)[ItVoc->first];
            nidkr = (Zeta[k].second)[ItVoc->first];
            for(size_t nid = 0; nid< nidkl+ nidkr; nid++){
                Data.push_back(ItVoc->first);                          // i dati sono scalati
            }

        }
    }

    file.write(reinterpret_cast<char*>(Data.data()),sizeof(POINT)*Data.size());
    file.close();

}
*/
/*
template<unsigned int DIM>
void CategoricalDocument<DIM>::ViewDocument(){

std::cout<<"Zeta.size() "<<Zeta.size()<<std::endl;


    unsigned int sum=0;
    std::cout<< "Alpha: "<<alpha<<std::endl;
	 std::cout<< "Nj: "<< Nj<<std::endl;

	if(Pi.empty())
	   std::cout<< "Pi is empty"<<std::endl;
	 std::cout<< "stampo Pi"<<std::endl;
	 for(auto i: Pi)
	    std::cout<< i<< "  ";
	std::cout<<std::endl;
	if(PiLeft.empty())
	   std::cout<< "PiLeft is empty"<<std::endl;
	std::cout<< "stampo PiLeft"<<std::endl;
	 for(auto i: PiLeft)
	    std::cout<< i<< "  ";
	std::cout<<std::endl;
	 if(PiRight.empty())
	   std::cout<< "PiRight is empty"<<std::endl;
	std::cout<< "stampo PiRight"<<std::endl;
	 for(auto i: PiRight)
	    std::cout<< i<< "  ";
	std::cout<<std::endl;
	if(LocalTable.empty())
	   std::cout<< "LocalTable is empty"<<std::endl;
	 std::cout<< "stampo LocalTable"<<std::endl;
	 for(auto i: LocalTable)
	    std::cout<< i<< "  ";
	std::cout<<std::endl;
	if(LocalTableLeft.empty())
	   std::cout<< "LocalTableLeft is empty"<<std::endl;
	std::cout<< "stampo LocalTableLeft"<<std::endl;
	for(auto i: LocalTableLeft)
	    std::cout<< i<< "  ";
	std::cout<<std::endl;
	if(LocalTableRight.empty())
	   std::cout<< "LocalTableRight is empty"<<std::endl;
	std::cout<< "stampo LocalTableRight"<<std::endl;
	for(auto i: LocalTableRight)
	    std::cout<< i<< "  ";
	std::cout<<std::endl;

    std::cout<< "Stampo il vocabolario:"<<std::endl;
    for(auto it= Vocabulary.cbegin(); it!= Vocabulary.cend() ; ++it){
       std::cout<< "key: " << it->first <<std::endl;
	   std::cout<< "value: "<<it->second<<std::endl;
	   }
	   
	std::cout<< "Stampo la Zeta:"<<std::endl;

	for(auto Zeta_it= Zeta.cbegin(); Zeta_it != Zeta.cend(); ++Zeta_it){
	   std::cout<< "Zeta key: "<< Zeta_it-> first<<std::endl;
	   std::cout<< "Zeta value: "<<std::endl;
	   auto left_map= Zeta_it->second.first;
	   auto right_map= Zeta_it->second.second;
	   std::cout<< "Left map: "<<std::endl;
	   for(auto left_it= left_map.cbegin(); left_it!= left_map.cend(); ++left_it){
	      std::cout<< "left map key, wordId: "<< left_it->first<<std::endl;
		  std::cout<< "left map value, nidjk: "<< left_it->second<<std::endl;
	   }
	   std::cout<< "Right map: "<<std::endl;
	    for(auto right_it= right_map.cbegin(); right_it!= right_map.cend(); ++right_it){
	      std::cout<< "right map key, wordId: "<< right_it->first<<std::endl;
		  std::cout<< "right map value, nidjk: "<< right_it->second<<std::endl;
	   }

	}
     sum=0;
	 std::cout<< "Stampo WordCount: "<<std::endl;
	 for(auto i: WordCount){
	    std::cout<< i<< "  ";
	    sum +=i;
    }
	 std::cout<<std::endl;
	 std::cout<< "Stampo WordCountLeft: "<<std::endl;
	 for(auto i: WordCountLeft)
	    std::cout<< i<< "  ";
	 std::cout<<std::endl;
	  std::cout<< "Stampo WordCountRight: "<<std::endl;
	 for(auto i: WordCountRight)
	    std::cout<< i<< "  ";
	 std::cout<<std::endl;
	 std::cout<<"SOMMA: "<<sum<<std::endl;

}
*/

#endif
