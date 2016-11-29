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
 * \brief  This file contains classes which manage the documents or, in general, groups of data.
 *  The generic class provides the common interface, whereas the derived and specialized classes are specific to the model.
 * \date February 2016
 */

 
/*! \brief Generic class for groups of data
 *
 *	Abstract class in which all methods are virtual.
 *  It contains methods that must be defined in all derived classes.
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */

 

template<typename Type, unsigned int DIM>
class GenericDocument{

	public:
        /*! \brief Updates group specific clusters' weights
		 *  \param _AllBeta - clusters' global weights
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdatePi (const vector<double>& , omprng&  ) = 0;

		/*! \brief 	Updates group specific weights for subclusters of cluster k
		 *  \param _BetaLeft - global weight for left subcluster of cluster k           
		 *  \param _BetaRight - global weight for right subcluster of cluster k
		 *  \param k - cluster id
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdatePiSub (const double , const double , const unsigned int , omprng& ) = 0;
		
				
	   	/*! \brief Updates group specific subclusters' weights for all cluster
		 *  \param _BetaLeft - global weights for left subclusters
		 *  \param _BetaRight - global weights for right subclusters
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdateAllPiSub(const vector<double> , const vector<double> , omprng& ) = 0;
		
		/*! \brief Updates table counts
		 *  \param _stirling - Stirling numbers
		 *  \param _Beta - clusters' global weights
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdateLocalTable (const vector<long double>& , const vector<double>& , omprng& ) = 0;
		
		/*! \brief Updates table counts for subclusters of cluster k
		 *  \param _stirling - Stirling numbers
		 *  \param _BetaLeft - global weight for left subcluster of cluster k
		 *  \param _BetaRight - global weight for right subcluster of cluster k
		 *  \param k - cluster id
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdateLocalTableSub_OneCluster(const vector<long double>& , const double _, const double , const unsigned int, omprng& ) = 0;
		
		/*! \brief Updates table counts for subclusters of all clusters
		 *  \param _stirling - Stirling numbers
		 *  \param _BetaLeft - global weight for left subclusters of all cluster
		 *  \param _BetaRight - global weight for right subclusters of all cluster 
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdateAllLocalTableSub (const vector<long double>& , const vector<double>& , const vector<double>& , omprng& ) = 0;
		
		 /*! \brief Updates cluster label of a datum, sampled with the Sampling method
		 *  \param _ThetaId - vector of a datum's weights in all clusters
		 *  \param _VettId - datum id
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdateZeta (const typename Type::THETA& , const unsigned int , omprng& )= 0;

		/*! \brief Updates cluster label of a datum, sampled with with the Sampling method
		 *  \param _ThetaId - vector of a datum's weights in all clusters
		 *  \param _ThetaIdLeft - vector of a datum's weights in all left subclusters
		 *  \param _ThetaIdRight - vector of a datum's weights in all right subclusters
		 *  \param _VettId - datum id
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdateZeta_and_Sub(const typename Type::THETA&, const typename Type::THETA& ,const typename Type::THETA& , const unsigned int , omprng& )= 0;

		/*! \brief Assigns the datum to the subclusters of cluster k, after sampling the subclusters' label with the Sampling method
		 *  \param _ThetaIdLeft - vector of a datum's weights in all left subclusters
		 *  \param _ThetaIdRight - vector of a datum's weights in all right subclusters
		 *  \param id - datum id
		 *  \param nidjk - number of times datum id in group j is assigned to cluster k
		 *  \param k - cluster id
		 *  \param Gen - parallel random number generator
		 */
		virtual void UpdateZetaSub(const typename Type::THETA& ,const typename Type::THETA&  ,const typename Type::Point , const unsigned int , const unsigned int , omprng& ) = 0;

		/*! \brief Method needed to remove cluster k during M-H moves
		 *  \param _k - cluster id
		 */
		virtual void UpdateZeta(const unsigned int ) = 0;

		/*! \brief Method needed to remove two clusters during M-H moves
		 *  \param _k1 - cluster id
		 *  \param _k2 - cluster id
		 */
		virtual void UpdateZeta(const unsigned int , const unsigned int ) = 0;

		/*! \brief Retrieve number of data in group j
		 *  \return Number of data in group j
		 */
		virtual unsigned int ViewNj() const = 0; 

	     /*! \brief Retrieve data id in group j
		 *  \param _VettId - filled with data id in group j
		 */
		virtual void ViewData( vector<typename Type::Point>& ) const = 0;      

		/*! \brief Retrieve counts necessary for updating subclusters' latent parameters
		 *  \param _k - cluster id
		 *  \param _counts4cleft - counts for left subcluster's latent parameter
		 *  \param _counts4cright - counts for right subcluster's latent parameter
		 */
		virtual void ViewCounts4c(const unsigned int , typename Type::STAT& , typename Type::STAT& )= 0; 

		/*! \brief Retrieve number of tables for cluster k \f$ m_{jk} \f$
		 *  \param _k - cluster id
		 *  \return \f$ m_{jk} \f$
		 */
		virtual unsigned int ViewNumTableID(const unsigned int) const = 0;

		/*! \brief Retrieve number of tables for left subcluster of cluster k \f$ m_{jkl} \f$
		 *  \param _k - cluster id
		 *  \return \f$ m_{jk} \f$
		 */
		virtual unsigned int ViewNumTableLeftID(const unsigned int ) const = 0;

		/*! \brief Retrieve number of tables for right subcluster of cluster k \f$ m_{jkr} \f$
		 *  \param _k - cluster id
		 *  \return \f$ m_{jkr} \f$
		 */
		virtual unsigned int ViewNumTableRightID(const unsigned int) const = 0;
		
		/*! \brief Retrieve number of data in cluster k \f$ n_{\cdot k}\f$
		 *  \param _k - cluster k
		 *  \return \f$ n_{\cdot k}\f$
		 */
		virtual unsigned int ViewDataCountID(const unsigned int ) const = 0;
		
		/*! \brief Retrieve number of data in left subcluster of cluster k \f$ n_{\cdot kl}\f$
		 *  \param _k - cluster id
		 *  \return \f$ n_{\cdot kl}\f$
		 */
		virtual unsigned int ViewDataCountLeftID(const unsigned int ) const = 0;

		/*! \brief Retrieve number of data in right subcluster of cluster k \f$ n_{\cdot kr}\f$
		 *  \param _k - cluster id
		 *  \return \f$ n_{\cdot kr}\f$
		 */
		virtual unsigned int ViewDataCountRightID(const unsigned int ) const = 0;
		
		/*! \brief Reset data counts in cluster k
		 *  \param k - cluster id
		 */
		virtual void ResetDataCountSub(const unsigned int) = 0;
		
		/*! \brief Retrieve a vector containing the number of data in all clusters
		 *  \param _WordCount - vector containing the number of data in all clusters
		 */
		virtual void ViewDataCount(vector<unsigned int>& ) const = 0;
		
		/*! \brief Retrieve a vector containing the number of data in all left subclusters
		 *  \param _WordCountLeft - vector containing the number of data in all left subclusters
		 */
		virtual void ViewDataCountLeft(vector<unsigned int>& ) const = 0;

		/*! \brief Retrieve a vector containing the number of data in all right subclusters
		 *  \param _WordCountRight - vector containing the number of data in all right subclusters
		 */
		virtual void ViewDataCountRight(vector<unsigned int>& ) const = 0;
		
		/*! \brief Retrieve counts and id of data in cluster k
		 *  \param _nidjk - structure that will contain counts and id of data in cluster k
		 *  \param _k - cluster id
		 */
		virtual void ViewIdCounts(vector<pair<typename Type::Point,unsigned int>>& ,const unsigned int ) = 0; 
		
		/*! \brief Retrieve cluster k
		 *  \param _k - cluster id
		 *  \param _Cluster - structure that will contain cluster k
		 */
		virtual void ViewCluster(const unsigned int, pair<unordered_map<typename Type::Point,unsigned int>,unordered_map<typename Type::Point,unsigned int>>& ) = 0;  
		
		/*! \brief Retrieve group specific weight of cluster k \f$ \pi_{jk} \f$
		 *  \param _k - cluster id
		 *  \return \f$ \pi_{jk} \f$
		 */
		virtual double ViewPiID(const unsigned int )const = 0;
		
		/*! \brief Retrieve group specific weight for the left subcluster of cluster k \f$ \pi_{jkl} \f$
		 *  \param _k - cluster id
		 *  \return \f$ \pi_{jkl} \f$
		 */
		virtual double ViewPiLeftID(const unsigned int ) const = 0;

		/*! \brief Retrieve group specific weight for the right subcluster of cluster k \f$ \pi_{jkr} \f$
		 *  \param _k - cluster id
		 *  \return \f$ \pi_{jkr} \f$
		 */
		virtual double ViewPiRightID(const unsigned int) const = 0;

		/*! \brief Retrieve a vector containing the group specific cluster weights
		 *  \param _pi - vector containing the group specific cluster weights
		 */
		virtual void ViewPi( vector<double>& ) const = 0;

		/*! \brief Retrieve a vector containing the group specific weights for left subclusters
		 *  \param _pi_left - vector containing the group specific weights for left subclusters
		 */
		virtual void ViewPiLeft( vector<double>& ) const = 0;
		
		/*! \brief Retrieve a vector containing the group specific weights for right subclusters
		 *  \param _pi_right - vector containing the group specific weights for right subclusters
		 */
		virtual void ViewPiRight( vector<double>& ) const = 0;

		/*! \brief Set the concentration parameter \f$ \alpha \f$ of the Dirichlet process governing the group
		 *  \param _alpha - \f$ \alpha \f$
		 */
		virtual void SetAlpha (const double ) = 0;

		/*! \brief Set the number of data in group j
		 *  \param _Nj - number of data in group j
		 */
		virtual void SetNj (const unsigned int ) = 0;
		
		/*! \brief Set the vector containing the group specific cluster weights
		 *  \param _pi - vector containing the group specific cluster weights
		 */
		virtual void SetPi(vector<double>& ) = 0;

	    /*! \brief Insert a new cluster
		 *  \param NewCluster - the new cluster
		 *  \param _Pi - group specific weight for the new cluster
		 *  \param _PiLeft - group specific left subcluster weight for the new cluster
		 *  \param _PiRight - group specific right subcluster weight for the new cluster
		 *  \param _WordCount - number of data in the new cluster
		 *  \param _WordCountLeft - number of data in the left subcluster of the new cluster
		 *  \param _WordCountRight - number of data in the right subcluster of the new cluster
		 *  \param _LocalTable - number of tables serving the new dish in restaurant j (CRF metaphor) \f$ m_{jk} \f$
		 *  \param _LocalTableLeft - number of tables serving the new left dish in restaurant j (CRF metaphor) \f$ m_{jkl} \f$
		 *  \param _LocalTableRight - number of tables serving the new right dish in restaurant j (CRF metaphor) \f$ m_{jkr} \f$
		 */
		virtual void InsertNewCluster(const pair < unordered_map < typename Type::Point, unsigned int>, unordered_map < typename Type::Point, unsigned int> >& , const double, const double, const double ,
	                   					   const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int )= 0;

		/*! \brief Remove clusters with id contained in the input vector
		 *  \param _k - vector containing the id of clusters to eliminate
		 */
		virtual void RemoveCluster(const vector<unsigned int>& ) = 0;
		
		/*! \brief Remove cluster k
		 *  \param _k - cluster id
		 */
		virtual void RemoveCluster( const unsigned int ) = 0;
	
		/*! \brief Remove two clusters
		 *  \param _k1 - cluster id
		 *  \param _k2 - cluster id
		 */
		virtual void RemoveCluster(unsigned int ,unsigned int ) = 0;

		/*! \brief Verifies if cluster k has an empty left subcluster
		 *  \param _k - cluster id
		 */
		virtual unsigned int CheckLeftSubcluster(const unsigned int ) = 0; 

		/*! \brief Verifies if cluster k has an empty right subcluster
		 *  \param _k - cluster id
		 */
		virtual unsigned int CheckRightSubcluster(const unsigned int ) = 0;  

		/*! \brief Retrieve labels assigned to data
		 *  \param Data - structure to store retrieved labels
		 */
		virtual void ViewLabel(vector<pair<typename Type::Point,unsigned int>>& ) = 0;    
	
		/*! \brief Acquires data
		 *  \param SSTR - contains data id and counts in the group
		 */   
		virtual void SetDataset(std::istringstream&) = 0;
	
        /*! \brief Allocates data in the Zeta container
		 *  \param _K - initial number of clusters
		 *  \param Gen - parallel random number generator
		 */	
		virtual unsigned int SortData(unsigned int, omprng&) = 0;


    private:

	   	/*! \brief Updates data counts in clusters
		 */
		virtual void UpdateDataCount() = 0;
	
		/*! \brief Samples cluster's or subcluster's label from categorical distribution
		 *  \param _temp_counts - vector containing count for datum id in all clusters
		 *  \param _Weights - weights for sampling labels
		 *  \param _nidj - number of times datum id appears in group j
		 *  \param Gen - parallel random number generator
		 */
		virtual void Sampling (std::vector<unsigned int>& , std::vector<double>& , unsigned int, omprng& ) = 0;   

};


/*! \brief Derived class for topic modeling, where data are categorical and the base measure is the Dirichlet distribution.
 * This class represents a document in the topic modeling problem.
 * It manages the words and is in charge of sampling the topic labels.
 * It samples the model's parameters specific to the document: \f$ \alpha, \pi_{j}, \bar(\pi)_{jl}, \bar(\pi)_{jr}, m_{j}, \bar{m}_{jl}, \bar{m}_{jl} \f$.
 * It keep track of words' counts in topics. 
 *
 * \authors{Debora Parisi and Stefania Perego}
 * \date Febbraio 2016
 */


template <unsigned int DIM=1>
class CategoricalDocument final: public GenericDocument <TypeCategorical<DIM>,DIM>{

   	public:
	
		 /*! \brief Statistics for updating hyperparameters of latent parameter's distribution
         */
		using STAT  = TypeCategorical<1>::STAT;
		
		/*! \brief Latent parameter: vector of distinct words' weights in the topic
         */   
		using THETA = TypeCategorical<1>::THETA;
		
		/*! \brief A datum
        */
		using POINT = TypeCategorical<1>::Point;
		
		/*! \brief Topic id
        */
		using ClusterID = unsigned int;
	
    private:

		/*! \brief Data container.
		 * Each topic has a left and right subtopic: in the subtopic map, the key is the datum and the mapped value is the number of times that datum appear in the document,
		 * in the subtopic of that topic.
         */
		 unordered_map < ClusterID , pair < unordered_map < POINT, unsigned int>, unordered_map < POINT, unsigned int> > >  Zeta;

	    /*! \brief Vocabulary of distinct words in the document.
		 *  The map key is the datum, the mapped value is the number of times the datum appears in the document.
		 */
		unordered_map <POINT, unsigned int> Vocabulary;
		
		/*! \brief Concentration parameter for the Dirichlet process ruling the document
		 */
		double alpha;
		
		/*! \brief Number of words contained in the document
		 */
		unsigned int Nj;
		
		/*! \brief Vector of group specific topics' weights
         */		
		vector<double> Pi;
		
		/*! \brief Vector of group specific weights for left subtopics
		*/	
		vector<double> PiLeft;
		
		/*! \brief Vector of group specific weights for right subtopics
         */	
		vector<double> PiRight;

        /*! \brief K-dimensional vector of counts: the k-th element represents the number of data in topic k (\f$ n_{jk} \f$)
         */		
		vector<unsigned int> WordCount;
		
		/*! \brief K-dimensional vector of counts: the k-th element represents the number of data in the left subtopic of topic k (\f$ n_{jkl} \f$)
         */	
		vector<unsigned int> WordCountLeft;
		
		/*! \brief K-dimensional vector of counts: the k-th element represents the number of data in the right subtopic of topic k (\f$ n_{jkr} \f$)
         */
		vector<unsigned int> WordCountRight;
		
		/*! \brief K-dimensional vector for tables: the k-th element represents the number of table in the restaurant serving dish k (\f$ m_{jk} \f$)
		 */
		vector<unsigned int> LocalTable;
		
		/*! \brief K-dimensional vector for tables: the k-th element represents the number of table in the restaurant serving dish k left (\f$ m_{jkl} \f$)
		 */
		vector<unsigned int> LocalTableLeft;
		
		/*! \brief K-dimensional vector for tables: the k-th element represents the number of table in the restaurant serving dish k right (\f$ m_{jkr} \f$)
		 */
		vector<unsigned int> LocalTableRight;

	public:
	
		/*! \brief Constructor
	     *  \param _alpha - concentration parameter of the Dirichlet process governing the document
		 */  
		CategoricalDocument(double _alpha): alpha(_alpha), Nj(0){};
		
		/*! \brief Default constructor
		 */
	    CategoricalDocument () = default;
		
		/*! \brief Destructor
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

        /*! \brief Update document specific topics' weights
		 *  \param _AllBeta - global topics' weights
		 *  \param Gen - parallel random number generator
		 */
		 
		void UpdatePi (const vector<double>& _AllBeta, omprng& Gen );
		
		/*! \brief 	Updates document specific weights for subtopics of topic k
		 *  \param _BetaLeft - global weight for left subtopic of topic k           
		 *  \param _BetaRight - global weight for right subtopic of topic k
		 *  \param k - topic id
		 *  \param Gen - parallel random number generator
		 */
		void UpdatePiSub (const double _BetaLeft, const double _BetaRight, const ClusterID k, omprng& Gen); 
		
	   	/*! \brief Updates document specific subtopics' weights for all topics
		 *  \param _BetaLeft - global weights for left subtopics
		 *  \param _BetaRight - global weights for right subtopics
		 *  \param Gen - parallel random number generator
		 */
	    void UpdateAllPiSub(const vector<double> _BetaLeft, const vector<double> _BetaRight, omprng& Gen); 

		/*! \brief Updates table counts
		 *  \param _stirling - Stirling numbers
		 *  \param _Beta - topics' global weights
		 *  \param Gen - parallel random number generator
		 */
    	void UpdateLocalTable(const vector<long double>& _stirling, const vector<double>& _Beta, omprng& Gen);
		
		/*! \brief Updates table counts for subtopics of topic k
		 *  \param _stirling - Stirling numbers
		 *  \param _BetaLeft - global weight for left subtopic of topic k
		 *  \param _BetaRight - global weight for right subtopic of topic k
		 *  \param k - topic id
		 *  \param Gen - parallel random number generator
		 */
		void UpdateLocalTableSub_OneCluster(const vector<long double>& _stirling, const double _BetaLeft, const double _BetaRight, const ClusterID k, omprng& Gen);
       
		/*! \brief Updates table counts for subtopic of all topics
		 *  \param _stirling - Stirling numbers
		 *  \param _BetaLeft - global weight for left subtopics of all topics
		 *  \param _BetaRight - global weight for right subtopics of all topics
		 *  \param Gen - parallel random number generator
		 */
     	void UpdateAllLocalTableSub(const vector<long double>& _stirling, const vector<double>& _BetaLeft, const vector<double>& _BetaRight, omprng& Gen);

		 /*! \brief Updates topic label of a datum, sampled with the Sampling method
		 *  \param _ThetaId - vector of a datum's weights in all topics
		 *  \param _VettId - datum id
		 *  \param Gen - parallel random number generator
		 */
	    void UpdateZeta(const THETA& _ThetaId, const POINT _VettId, omprng& Gen);
    
		/*! \brief Updates topic label of a datum, sampled with with the Sampling method
		 *  \param _ThetaId - vector of a datum's weights in all topics
		 *  \param _ThetaIdLeft - vector of a datum's weights in all left subtopics
		 *  \param _ThetaIdRight - vector of a datum's weights in all right subtopics
		 *  \param _VettId - datum id
		 *  \param Gen - parallel random number generator
		 */
		void UpdateZeta_and_Sub(const THETA& _ThetaId, const THETA& _ThetaIdLeft,const THETA& _ThetaIdRight, const unsigned int _VettId, omprng& Gen);
	    
		/*! \brief Assigns the datum to the subtopics of topic k, after sampling the subtopics' label with the Sampling method
		 *  \param _ThetaIdLeft - vector of a datum's weights in all left subtopics
		 *  \param _ThetaIdRight - vector of a datum's weights in all right subtopics
		 *  \param id - datum id
		 *  \param nidjk - number of times datum id in document j is assigned to topic k
		 *  \param k - topic id
		 *  \param Gen - parallel random number generator
		 */
		void UpdateZetaSub(const THETA& _ThetaIdLeft,const THETA& _ThetaIdRight ,const POINT id, const unsigned int nidjk, const ClusterID k, omprng& Gen);

		/*! \brief Method needed to remove topic k during M-H moves
		 *  \param _k - topic id
		 */
		void UpdateZeta(const ClusterID _k);
		
		/*! \brief Method needed to remove two topics during M-H moves
		 *  \param _k1 - topic id
		 *  \param _k2 - topic id
		 */
		void UpdateZeta(const ClusterID _k1, const ClusterID _k2);

		/*! \brief Retrieve number of data in document j
		 *  \return Number of data in document j
		 */
		unsigned int ViewNj() const ;
		
	     /*! \brief Retrieve data id in document j
		 *  \param _VettId - filled with data id in document j
		 */
		void ViewData( vector<POINT>& _VettId) const;
		
		/*! \brief Retrieve counts necessary for updating subtopics' latent parameters
		 *  \param _k - topic id
		 *  \param _counts4cleft - counts for left subtopic's latent parameter
		 *  \param _counts4cright - counts for right subtopic's latent parameter
		 */
		void ViewCounts4c(ClusterID _k, STAT& _counts4cleft, STAT& _counts4cright);


		/*! \brief Retrieve number of tables for topic k \f$ m_{jk} \f$
		 *  \param _k - topic id
		 *  \return \f$ m_{jk} \f$
		 */
		unsigned int ViewNumTableID(const ClusterID _k) const;
		
		/*! \brief Retrieve number of tables for left subtopic of topic k \f$ m_{jkl} \f$
		 *  \param _k - topic id
		 *  \return \f$ m_{jk} \f$
		 */
        unsigned int ViewNumTableLeftID(const ClusterID _k) const;
	
		/*! \brief Retrieve number of tables for right subtopic of topic k \f$ m_{jkr} \f$
		 *  \param _k - topic id
		 *  \return \f$ m_{jkr} \f$
		 */
		unsigned int ViewNumTableRightID(const ClusterID _k) const;

		/*! \brief Retrieve number of words in topic k \f$ n_{\cdot k}\f$
		 *  \param _k - topic k
		 *  \return \f$ n_{\cdot k}\f$
		 */
	    unsigned int ViewDataCountID(const ClusterID _k) const;
		
		/*! \brief Retrieve number of topic in left subtopic of topic k \f$ n_{\cdot kl}\f$
		 *  \param _k - topic id
		 *  \return \f$ n_{\cdot kl}\f$
		 */
		unsigned int ViewDataCountLeftID(const ClusterID _k)const;
		
		/*! \brief Retrieve number of topic in right subtopic of topic k \f$ n_{\cdot kr}\f$
		 *  \param _k - topic id
		 *  \return \f$ n_{\cdot kr}\f$
		 */
		unsigned int ViewDataCountRightID(const ClusterID _k)const;
		
		/*! \brief Reset wourds counts in topic k
		 *  \param k - topic id
		 */
		void ResetDataCountSub(const ClusterID k);

		/*! \brief Retrieve a vector containing the number of words in all topics
		 *  \param _WordCount - vector containing the number of words in all 
		 */
		void ViewDataCount(vector<unsigned int>& _WordCount) const;
		
		
		/*! \brief Retrieve a vector containing the number of words in all left subtopics
		 *  \param _WordCountLeft - vector containing the number of words in all left subtopics
		 */
		void ViewDataCountLeft(vector<unsigned int>& _WordCountLeft) const;
		
		
		/*! \brief Retrieve a vector containing the number of words in all right subtopics
		 *  \param _WordCountRight - vector containing the number of words in all right subtopics
		 */
		void ViewDataCountRight(vector<unsigned int>& _WordCountRight) const;
		
		/*! \brief Retrieve counts and id of words in topic k
		 *  \param _nidjk - structure that will contain counts and id of words in topic k
		 *  \param _k - topic id
		 */
		void ViewIdCounts(vector<pair<POINT,unsigned int>>& _nidjk,const ClusterID _k);

		/*! \brief Retrieve topic k
		 *  \param _k - topic id
		 *  \param _Cluster - structure that will contain topic k
		 */ 
		void ViewCluster(const ClusterID _k, pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>>& _Cluster) ;        

		/*! \brief Retrieve document specific weight of topic k \f$ \pi_{jk} \f$
		 *  \param _k - topic id
		 *  \return \f$ \pi_{jk} \f$
		 */
		double ViewPiID(const ClusterID _k) const;
		
		/*! \brief Retrieve document specific weight for the left subtopic of topic k \f$ \pi_{jkl} \f$
		 *  \param _k - topic id
		 *  \return \f$ \pi_{jkl} \f$
		 */
		double ViewPiLeftID(const ClusterID _k) const;
		
		/*! \brief Retrieve document specific weight for the right subtopic of topic k \f$ \pi_{jkr} \f$
		 *  \param _k - topic id
		 *  \return \f$ \pi_{jkr} \f$
		 */
		double ViewPiRightID(const ClusterID _k) const;		
		
		/*! \brief Retrieve a vector containing the document specific topic weights
		 *  \param _pi - vector containing the document specific topic weights
		 */
		void ViewPi( vector<double>& _pi) const;
		
		/*! \brief Retrieve a vector containing the document specific weights for left subtopics
		 *  \param _pi_left - vector containing the document specific weights for left subtopics
		 */
		void ViewPiLeft( vector<double>& _pi) const;
		
		/*! \brief Retrieve a vector containing the document specific weights for right subtopics
		 *  \param _pi_right - vector containing the document specific weights for right subtopics
		 */
		void ViewPiRight( vector<double>& _pi) const;
		
		/*! \brief Set the concentration parameter \f$ \alpha \f$ of the Dirichlet process governing the document
		 *  \param _alpha - \f$ \alpha \f$
		 */
		void SetAlpha (const double _alpha);
		
		/*! \brief Set the number of words in document j
		 *  \param _Nj - number of words in document j
		 */
		void SetNj (const unsigned int _Nj);
		
		/*! \brief Set the vector containing the document specific topics weights
		 *  \param _pi - vector containing the document specific topics weights
		 */
		void SetPi(vector<double>& _Pi);

	    /*! \brief Insert a new topic
		 *  \param NewCluster - the new topic
		 *  \param _Pi - document specific weight for the new topic
		 *  \param _PiLeft - document specific left subtopic weight for the new topic
		 *  \param _PiRight - document specific right subtopic weight for the new topic
		 *  \param _WordCount - number of words in the new topic
		 *  \param _WordCountLeft - number of words in the left subtopic of the new topic
		 *  \param _WordCountRight - number of words in the right subtopic of the new topic
		 *  \param _LocalTable - number of tables serving the new dish in restaurant j (CRF metaphor) \f$ m_{jk} \f$
		 *  \param _LocalTableLeft - number of tables serving the new left dish in restaurant j (CRF metaphor) \f$ m_{jkl} \f$
		 *  \param _LocalTableRight - number of tables serving the new right dish in restaurant j (CRF metaphor) \f$ m_{jkr} \f$
		 */
		void InsertNewCluster(const pair < unordered_map < POINT, unsigned int>, unordered_map < POINT, unsigned int> >& NewCluster, const double _Pi, const double _PiLeft, 
	                      const double _PiRight, const unsigned int _WordCount, const unsigned int _WordCountLeft, const unsigned int _WordCountRight, const unsigned int _LocalTable,
						  const unsigned int _LocalTableLeft, const unsigned int _LocalTableRight);
						  
		/*! \brief Remove topics with id contained in the input vector
		 *  \param _k - vector containing the id of topics to eliminate
		 */
		void RemoveCluster(const vector<ClusterID>& _k);
		
		/*! \brief Remove topic k
		 *  \param _k - topic id
		 */
		void RemoveCluster(const  ClusterID _k);
		
		/*! \brief Remove two topics
		 *  \param _k1 - topic id
		 *  \param _k2 - topic id
		 */
		void RemoveCluster(ClusterID _k1, ClusterID _k2);

		/*! \brief Verifies if topic k has an empty left subtopic
		 *  \param _k - topic id
		 */
		unsigned int CheckLeftSubcluster(const ClusterID _k);
		
		/*! \brief Verifies if topic k has an empty right subtopic
		 *  \param _k - topic id
		 */
		unsigned int CheckRightSubcluster(const ClusterID _k);

		/*! \brief Retrieve labels assigned to words
		 *  \param Data - structure to store retrieved labels
		 */
		void ViewLabel(vector<pair<POINT,ClusterID>>& Data);
	
		/*! \brief Acquires data
		 *  \param SSTR - contains words id and counts in the document
		 */
		void SetDataset(std::istringstream& SSTR);
		
        /*! \brief Allocates words in the Zeta container
		 *  \param _K - initial number of topics
		 *  \param Gen - parallel random number generator
		 */
		unsigned int SortData(unsigned int _K, omprng& Gen);


    private:
		
	   	/*! \brief Updates word counts in topics
		 */
		void UpdateDataCount();
		
		/*! \brief Samples topic's or subtopic's label from categorical distribution
		 *  \param _temp_counts - vector containing count for word id in all topics
		 *  \param _Weights - weights for sampling labels
		 *  \param _nidj - number of times word id appears in document j
		 *  \param Gen - parallel random number generator
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
        std::cerr<<"Error in Update Zeta and Sub"<<std::endl;
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


#endif
