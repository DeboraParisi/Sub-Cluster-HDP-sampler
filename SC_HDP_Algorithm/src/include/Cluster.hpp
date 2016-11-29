#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <vector>
#include <iostream>
#include "Type.hpp"

using std::vector;

/*! \file Cluster.hpp
 *
 * \brief Data's structures which manage the cluster. This structure depend on the model.
 * These classes define how to manage all parameter which are involve in the definition of cluster.
 * In these classes there aren't any methods which sample some variables that describe the cluster. There are only methods that read, write and 
 * keep in memory information of cluster and his sub-clusters.
 * \date February 2016
 */


/*! \brief Generic Model of Cluster
 *
 *  Abstract class where all methods are virtual and they are null.
 *	All classes that inherit form Cluster Generic are used to extract, to memorize and to set all data which describe a cluster and its subclusters.
 * 	This classes have the role only to manage the clusters, in this class there anern't any sample function.
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */

template <template <unsigned int> class ClassType, unsigned int DIM>
class GenericCluster
{
public:


    /*! \brief  It sets the latent parameter of cluster
     *  \param   _Theta -  latent parameter in input. Its type is THETA
     */
    virtual void SetTheta(typename ClassType<DIM>::THETA&) = 0;

    /*! \brief  It sets the left sub-cluster latent parameter
     *  \param  _ThetaLeft - latent parameter in input. Its type is THETA
     */
    virtual void SetThetaLeft(typename ClassType<DIM>::THETA&) = 0;

    /*! \brief  It sets the right sub-cluster latent parameter
     *  \param  _ThetaRight - latent parameter in input. Its type is THETA
     */
    virtual void SetThetaRight(typename ClassType<DIM>::THETA&) = 0;

    /*! \brief  It sets the weigth of cluster
     *  \param  _Beta - Weigth in input
     */
    virtual void SetBeta(double) = 0;

    /*! \brief	It sets the weigth of left subcluster
     *  \param _BetaLeft - Weigth in input
     */
    virtual void SetBetaLeft(double) = 0;

    /*! \brief	It sets the weigth of right subcluster
     *  \param _BetaRight - Weigth in input
     */
    virtual void SetBetaRight(double) = 0;

    /*! \brief It sets the number of global tables which characterize the cluster
        \param NrTable - number of tables which characterize the cluster
	 */
    virtual void SetGlobalTable (unsigned int) = 0;

    /*! \brief It sets the number of global tables which characterize the left sub-cluster 
        \param NrTableLeft - number of tables which characterize the left sub-cluster 
	 */
    virtual void SetGlobalTableLeft (unsigned int) = 0;

    /*! \brief It sets the number of global tables which characterize the right sub-cluster 
        \param NrTableRight - number of tables which characterize the right sub-cluster 
	 */
    virtual void SetGlobalTableRight (unsigned int) = 0;

	/*! \brief  It sets the statistics, in other words the hyperparameter of cluster's latent paramenter
     *  \param  _c - cluster's statistics
	 */
	virtual void SetStatistics(typename ClassType<DIM>::STAT&) = 0;

    /*! \brief  It sets the statistics, in other words the hyperparameter of left sub-cluster's latent paramenter
     *  \param  _cLeft - left sub-cluster's statistics
	 */
    virtual void SetStatisticsLeft(typename ClassType<DIM>::STAT&) = 0;

    /*! \brief  It sets the statistics, in other words the hyperparameter of right sub-cluster's latent paramenter
     *  \param _cRight - right sub-cluster's statistics
	 */
    virtual void SetStatisticsRight(typename ClassType<DIM>::STAT&) = 0;

    /*! \brief It extracts information about cluster's latent parameter
        \param _Theta - It takes in input an object, which its type is THETA. In this object It saves the value of cluster's latent paramenter
     */
    virtual void ViewTheta (typename ClassType<DIM>::THETA&) const = 0;

    /*! \brief It extracts information about left sub-cluster's latent parameter
        \param _ThetaLeft - It takes in input an object, which its type is THETA. In this object It saves the value of left sub-cluster's latent paramenter
     */
    virtual void ViewThetaLeft(typename ClassType<DIM>::THETA&) const = 0;

    /*! \brief It extracts information about right sub-cluster's latent parameter
        \param _ThetaRight - It takes in input an object, which its type is THETA. In this object It saves the value of right sub-cluster's latent paramenter
     */
    virtual void ViewThetaRight(typename ClassType<DIM>::THETA&) const = 0;

    /*! \briefIt It extracts information about cluster's global weight
     *  \return cluster's global weight
     */
    virtual double ViewBeta() const = 0;

    /*! \briefIt It extracts information about left sub-cluster's global weight
     *  \return left sub-cluster's global weight
     */
    virtual double ViewBetaLeft() const = 0;

    /*! \briefIt It extracts information about right sub-cluster's global weight
     *  \return right sub-cluster's global weight
     */
    virtual double ViewBetaRight() const = 0;

    /*! \brief It extracts infromation about statistics of cluster's latent paramenter
     *  \param _c - It takes in input an object, which its type is STAT. In this object It saves the statistics of cluster's latent paramenter
     */
    virtual void ViewStatistics( typename ClassType<DIM>::STAT&) const = 0;

    /*! \brief It extracts infromation about statistics of left sub-cluster's latent paramenter 
     *  \param _cLeft - It takes in input an object, which its type is STAT. In this object It saves the statistics of left sub-cluster's latent paramenter
     */
    virtual void ViewStatisticsLeft( typename ClassType<DIM>::STAT&) const = 0;

    /*! \brief It extracts infromation about statistics of right sub-cluster's latent paramenter 
     *  \param _cRight - It takes in input an object, which its type is STAT. In this object It saves the statistics of right sub-cluster's latent paramenter
     */
    virtual void ViewStatisticsRight( typename ClassType<DIM>::STAT&) const = 0;

    /*! \brief It extracts the number of cluster's global table 
     *  \return Number of global table
     */
    virtual unsigned int ViewGlobalTable () const = 0;

    /*! \brief It extracts the number of left sub-cluster's global table 
     *  \return Number of global table in left sub-cluster
     */
    virtual unsigned int ViewGlobalTableLeft () const = 0;

     /*! \brief It extracts the number of right sub-cluster's global table 
     *  \return Number of global table in right sub-cluster
     */
    virtual unsigned int ViewGlobalTableRight () const = 0;

    /*! \brief  It resests the cluster's statistics
     *  \param  W - dimenstion of statistics with to update hyperparamenters of cluster's latent paramenter 
	 */
	virtual void ResetStatistics(unsigned int) = 0;

    /*! \brief  It resests the left sub-cluster's statistics
     *  \param  W - dimenstion of statistics with to update hyperparamenters of left sub-cluster's latent paramenter 
	 */
	virtual void ResetStatisticsLeft(unsigned int) = 0;

	/*! \brief  It resests the right sub-cluster's statistics
     *  \param  W - dimenstion of statistics with to update hyperparamenters of right sub-cluster's latent paramenter 
	 */
	virtual void ResetStatisticsRight(unsigned int) = 0;

	/*! \brief  It updates the statistics, in other words hyperparamenters of cluster and sub-clusters 
     *  \param  counts4cleft - statistics to update the hyperparameters of left sub-subcluster's latent parameters
     *  \param  counts4cright - statistics to update the hyperparameters of right sub-subcluster's latent parameters
	 */
	virtual void UpdateStatistics(typename ClassType<DIM>::STAT&, typename ClassType<DIM>::STAT&) = 0;

    /*! \brief Check if cluster is empty
     *  \return TRUE if cluster is empty otherwise  FALSE 
     */
    virtual bool IsEmpty() const = 0;

};

/*! \brief Management of cluster's and subclusters' informations with Categorical Likelihood
 *
 *  This class memorizes and extracts information of cluster and sub-clusters about global weight, latent parameters, hyperparameters' update.
 * 	In this case that the likelihood is categorical, data could be repeted, the latent parameters are the weight of distinct elements and the statistics are conts of elements in the clusters.
 *	The latent paramenters are the mixutre's parameters.
 *	\authors{Debora Parisi and Stefania Perego}
 *  \date Febbrario 2016
 */
class CategoricalCluster final: GenericCluster<TypeCategorical,1>
{
    public:

        //DEFINITION OF PUBLIC TYPE

        /*! \brief Latent Paramenter, vector of distinct element's weigths which are contained into the cluster
         */
        using THETA = TypeCategorical<1>::THETA;

        /*! \brief Type of single data, repeaterd data.
         */
        using Point = TypeCategorical<1>::Point;

        /*! \brief Statistics' vector. In this container there are the update of latent parameters' hyperparameters of clusters or sub-clusters. Number of data which are contained into cluster and sub-clusters
         */
        using STAT  = TypeCategorical<1>::STAT;

    private:

        // DEFINITION OF PRIVATE TYPE

        /*! \brief Cluster's global weigth
         */
        double Beta;

        /*! \brief Left sub-cluster's global weigth
         */
        double BetaLeft;

         /*! \brief Right sub-cluster's global weigth
         */
        double BetaRight;

        /*! \brief Cluster's latent parameter: weigth of distnict element which are conteined into the cluster
         */
        THETA Theta;

        /*! \brief Left sub-cluster's latent parameter: weigth of distnict element which are conteined into the left sub-cluster
         */
        THETA ThetaLeft;

        /*! \brief Right sub-cluster's latent parameter: weigth of distnict element which are conteined into the right sub-cluster
         */
        THETA ThetaRight;

        /*! \brief Statisitcs to update clustr's latent paramenter: counts of element which are conteined into the cluster
         */
        STAT c;

        /*! \brief Statisitcs to update left sub-clustr's latent paramenter: counts of element which are conteined into the left sub-cluster
         */
        STAT cLeft;

        /*! \brief Statisitcs to update right sub-clustr's latent paramenter: counts of element which are conteined into the right sub-cluster
         */
        STAT cRight;

        /*! \brief Number of global table assings to the cluster
         */
        unsigned int NrTable;

         /*! \brief Number of global table assings to the left sub-cluster
         */
        unsigned int NrTableLeft;

         /*! \brief Number of global table assings to the right sub-cluster
         */
        unsigned int NrTableRight;

    public:

        //COSTRUCTOR and DISTRUCTOR

        /*! \brief Default costructor 
         */
        CategoricalCluster();

        /*! \brief  Default distructor 
         */
        ~CategoricalCluster()=default;

        /*! \brief Costructor which required in input all informations about cluster and sub-clusters
         *  \param _Beta - Cluster's global weigth
         *  \param _BetaLeft - Left sub-luster's global weigth
         *  \param _BetaRight - Right sub-luster's global weigth
         *  \param _Theta - Cluster's latent paramenter: weigth of distnict element which are conteined into the cluster
         *  \param _ThetaLeft - Left sub-cluster's latent parameter: weigth of distnict element which are conteined into the left sub-cluster
         *  \param _ThetaRight - Right sub-cluster's latent parameter: weigth of distnict element which are conteined into the right sub-cluster
         *  \param _c - Statisitcs to update clustr's latent paramenter: counts of element which are conteined into the cluster
         *  \param _cLeft -Statisitcs to update left sub-clustr's latent paramenter: counts of element which are conteined into the left sub-cluster
         *  \param _cRight - Statisitcs to update right sub-clustr's latent paramenter: counts of element which are conteined into the right sub-cluster
         *  \param -NrTable - Number of global table assings to the cluster
         *  \param -NrTableLeft - Number of global table assings to the left sub-cluster
         *  \param -NrTableRight - Number of global table assings to the right sub-cluster
         */
        CategoricalCluster(double _Beta, double _BetaLeft, double _BetaRight, THETA& _Theta ,THETA& _ThetaLeft, THETA& _ThetaRight, STAT& _c, STAT& _cLeft, STAT& _cRight, unsigned int _NrTable, unsigned int _NrTableLeft, unsigned int _NrTableRight);

       //METODI PER FISSARE
        /*! \brief  Fix the cluster's latent paramenter 
         *  \param   _Theta - cluster's latent paramenter. Its type is THETA
         */
        void SetTheta(THETA& _Theta);

        /*! \brief  Fix the left sub-cluster's latent paramenter 
         *  \param   _Theta - left sub-cluster's latent paramenter. Its type is THETA
         */
        void SetThetaLeft(THETA& _ThetaLeft);

        /*! \brief  Fix the right sub-cluster's latent paramenter 
         *  \param   _Theta - Right sub-cluster's latent paramenter. Its type is THETA
         */
        void SetThetaRight(THETA& _ThetaRight);

        /*! \brief  Fix cluster's global weigth
         *  \param  _Beta - Weigth in input 
         */
        void SetBeta(double _Beta);

        /*! \brief  Fix left sub-cluster's global weigth
         *  \param  _Beta - Weigth in input 
         */
        void SetBetaLeft(double _BetaLeft);

        /*! \brief  Fix rigth sub-cluster's global weigth
         *  \param  _Beta - Weigth in input 
         */
        void SetBetaRight(double _BetaRight);

        /*! \brief Fix the number of global table which characterizes the cluster 
         *  \param NrTable - number of global table which characterizes the cluster 
         */
        void SetGlobalTable(unsigned int _NrTable);

       /*! \brief Fix the number of global table which characterizes the left sub-cluster 
         *  \param NrTable - number of global table which characterizes the left sub-cluster 
         */
        void SetGlobalTableLeft(unsigned int _NrTableLeft);

        /*! \brief Fix the number of global table which characterizes the right sub-cluster 
         *  \param NrTable - number of global table which characterizes the right sub-cluster 
         */
        void SetGlobalTableRight(unsigned int _NrTableRight);

        /*! \brief  Set the statistics. Fissa le statistiche, ovvero gli iperparametri dei parametri latenti del cluster
         *  \param  _c - statistiche del cluster, conteggi degli elementi finiti nel cluster
         */
        void SetStatistics(STAT& _c);

        /*! \brief  Fissa le statistiche, ovvero gli iperparametri dei parametri latenti del sub-cluster sinistro
         *  \param  _cLeft - statistiche del sub-cluster sinistro, conteggi degli elementi finiti nel sub-cluster sinistro
         */
        void SetStatisticsLeft(STAT& _cLeft);

        /*! \brief  Fissa le statistiche, ovvero gli iperparametri dei parametri latenti del sub-cluster destro
         *  \param  _cLeft - statistiche del sub-cluster destro, conteggi degli elementi finiti nel sub-cluster destro
         */
        void SetStatisticsRight(STAT& _cRight);

        //METODI PER VISUALIZZARE

        /*! \brief Estrae il parametro latente del cluster
         *  \param _Theta - oggetto di tipo THETA in cui viene memorizzato il parametro latente del cluster, peso degli elementi distinti nel cluster
         */
        void ViewTheta (THETA& _Theta) const;

        /*! \brief Estrae il parametro latente del sub-cluster sinistro
         *  \param _ThetaLeft - oggetto di tipo THETA in cui viene memorizzato il parametro latente del cluster, peso degli elementi distinti
                                nel sub-cluster sinistro
         */
        void ViewThetaLeft(THETA& _ThetaLeft) const;

        /*! \brief Estrae il parametro latente del sub-cluster destro
         *  \param _ThetaRight - oggetto di tipo THETA in cui viene memorizzato il parametro latente del cluster, peso degli elementi distinti
                                nel sub-cluster destro
         */
        void ViewThetaRight(THETA& _ThetaRight) const;

        /*! \brief Estrare l'i-esimo elemento del parametro latente del cluster
         *  \param Elemento da estrarre
         *  \return peso dell'elemento nella posizione indicata del cluster
         */
        double ViewThetaId (unsigned int _id) const;

        /*! \brief Estrare l'i-esimo elemento del parametro latente del sub-cluster sinistro
         *  \param Elemento da estrarre
         *  \return peso dell'elemento nella posizione indicata del sub-cluster sinistro
         */
        double ViewThetaLeftId (unsigned int _id) const;

        /*! \brief Estrare l'i-esimo elemento del parametro latente del sub-cluster destro
         *  \param Elemento da estrarre
         *  \return peso dell'elemento nella posizione indicata del subl-cluster destro
         */
        double ViewThetaRightId (unsigned int _id) const;

        /*! \brief Estra il peso globale del cluster
         *  \return Peso del cluster
         */
        double ViewBeta() const;

        /*! \brief Estra il peso globale del sub-cluster sinistro
         *  \return Peso del sub-cluster sinistro
         */
        double ViewBetaLeft() const;

        /*! \brief Estra il peso globale del sub-cluster destro
         *  \return Peso del sub-cluster destro
         */
        double ViewBetaRight() const;

        /*! \brief Estra il numero globale dei tavoli nel cluster
         *  \return Numero di tavoli nel cluster
         */
        unsigned int ViewGlobalTable () const;

        /*! \brief Estra il numero globale dei tavoli nel sub-cluster sinistro
         *  \return Numero di tavoli nel sub-cluster sinistro
         */
        unsigned int ViewGlobalTableLeft () const;

        /*! \brief Estra il numero globale dei tavoli nel sub-cluster destro
         *  \return Numero di tavoli nel sub-cluster destrp
         */
        unsigned int ViewGlobalTableRight () const;

        /*! \brief  It sets the statistics, in other words the hyperparameter of cluster's latent paramenter
     	 *  \param  _c - Object's type is STA. 
	 	 */
        void ViewStatistics( STAT& _c) const;

         /*! \brief  It sets the statistics, in other words the hyperparameter of left sub-cluster's latent paramenter
     	 *  \param  _c - Object's type is STA. 
	 	 */
        void ViewStatisticsLeft( STAT& _cLeft) const;

         /*! \brief  It sets the statistics, in other words the hyperparameter of right sub-cluster's latent paramenter
     	 *  \param  _c - Object's type is STA. 
	 	 */
        void ViewStatisticsRight( STAT& _cRight) const;

        /*! \brief  It resests the cluster's statistics
    	 *  \param  W - dimenstion of statistics with to update hyperparamenters of cluster's latent paramenter 
		 */
        void ResetStatistics(unsigned int W);

        /*! \brief  It resests the left sub-cluster's statistics
    	 *  \param  W - dimenstion of statistics with to update hyperparamenters of left sub-cluster's latent paramenter 
		 */
        void ResetStatisticsLeft(unsigned int W);

        /*! \brief  It resests the right sub-cluster's statistics
    	 *  \param  W - dimenstion of statistics with to update hyperparamenters of right sub-cluster's latent paramenter 
		 */
        void ResetStatisticsRight(unsigned int W);

        /*! \brief  Update of statistics. This method update the hyperparameter of cluster's and sub-clusters' latent paramenter.
         *  \param  counts4cleft - statistics to update the hyperparameters of left sub-subcluster's latent parameters
     	 *  \param  counts4cright - statistics to update the hyperparameters of right sub-subcluster's latent parameters
	 	 */
        void UpdateStatistics(STAT& counts4cleft, STAT& counts4right);

        /*! \brief Check if cluster is empty.
         *  \return TRUE if cluster is empty otherwise FALSE 
         */
        bool IsEmpty() const;

};

////////////////////////////////////////
// DEFINIZIONI DI CLUSTERCATEGORICAL ///
////////////////////////////////////////

//Costruttore di default
CategoricalCluster::CategoricalCluster():
Beta(1.0/3), BetaLeft(0.5), BetaRight(0.5), NrTable(0), NrTableLeft(0), NrTableRight(0){}

//Costruttore che richiede tutte le informazioni del cluster e subcluster
CategoricalCluster::CategoricalCluster(double _Beta, double _BetaLeft, double _BetaRight, THETA& _Theta ,THETA& _ThetaLeft, THETA& _ThetaRight, STAT& _c, STAT& _cLeft, STAT& _cRight, unsigned int _NrTable, unsigned int _NrTableLeft, unsigned int _NrTableRight):
Beta(_Beta), BetaLeft(_BetaLeft), BetaRight(_BetaRight), NrTable(_NrTable), NrTableLeft(_NrTableLeft), NrTableRight(_NrTableRight) {

	Theta= _Theta;
	ThetaLeft= _ThetaLeft;
	ThetaRight= _ThetaRight;
	c= _c;
	cLeft= _cLeft;
	cRight= _cRight;

}

//Fissa il parametro latente del cluster
void CategoricalCluster::SetTheta(THETA& _Theta){
    Theta=_Theta;
}

//Fissa il parametro latente del sub-cluster sinistro
void CategoricalCluster::SetThetaLeft(THETA& _ThetaLeft){
    ThetaLeft=_ThetaLeft;
}

//Fissa il parametro latente del sub-cluster destro
void CategoricalCluster::SetThetaRight(THETA& _ThetaRight){
    ThetaRight=_ThetaRight;
}

//Fissa il peso globale del cluster
void CategoricalCluster::SetBeta(double _Beta){
    Beta=_Beta;
}

//Fissa il peso globale del sub-cluster sinistro
void CategoricalCluster::SetBetaLeft(double _BetaLeft){
    BetaLeft=_BetaLeft;
}

//Fissa il peso globale del sub-cluster destro
void CategoricalCluster::SetBetaRight(double _BetaRight){
    BetaRight=_BetaRight;
}

//Fissa la statistica del cluster
void CategoricalCluster::SetStatistics(STAT& _c){
    c=_c;
}

//Fissa la statistica del sub-cluster sinistro
void CategoricalCluster::SetStatisticsLeft(STAT& _cLeft){
    cLeft=_cLeft;
}

//Fissa la statistica del sub-cluster destro
void CategoricalCluster::SetStatisticsRight(STAT& _cRight){
    cRight=_cRight;
}

//Fissa il numero di tavoli globali del cluster
void CategoricalCluster::SetGlobalTable(unsigned int _NrTable){
    NrTable=_NrTable;
}

//Fissa il numero di tavoli globali del sub-cluster sinistro
void CategoricalCluster::SetGlobalTableLeft(unsigned int _NrTableLeft){
    NrTableLeft=_NrTableLeft;
}

//Fissa il numero i tavoli globali del sub-cluster destro
void CategoricalCluster::SetGlobalTableRight(unsigned int _NrTableRight){
    NrTableRight=_NrTableRight;
}

//Estrae i parametri latenti del cluster
void CategoricalCluster::ViewTheta(THETA& _Theta)const{
    _Theta=Theta;
}

//Estrae i parametri latenti del sub-cluster sinistro
void CategoricalCluster::ViewThetaLeft(THETA& _ThetaLeft) const{
    _ThetaLeft=ThetaLeft;
}

//Estrae i parametri lantenti del sub-cluster destro
void CategoricalCluster::ViewThetaRight(THETA& _ThetaRight) const{
    _ThetaRight=ThetaRight;
}

//Estrare l'i-esimo elemento del parametro latente del cluster
double CategoricalCluster::ViewThetaId (unsigned int _id) const{

    if(_id > Theta.size()){
        std::cerr<<"Errore in CategoricalCluster.ViewThetaId"<<std::endl;
        std::cerr<<"Id eccede rispetto gli elemtni contenuti"<<std::endl;
        exit(1);
    }
    return Theta[_id];
}

//Estrare l'i-esimo elemento del parametro latente del sub-cluster sinistro
double CategoricalCluster::ViewThetaLeftId (unsigned int _id) const{

    if(_id > ThetaLeft.size()){
        std::cerr<<"Errore in CategoricalCluster.ViewThetaLeftId"<<std::endl;
        std::cerr<<"Id eccede rispetto gli elemtni contenuti"<<std::endl;
        exit(1);
    }
    return ThetaLeft[_id];
}

//Estrare l'i-esimo elemento del parametro latente del sub-cluster destro
double CategoricalCluster::ViewThetaRightId (unsigned int _id) const{

    if(_id > ThetaRight.size()){
        std::cerr<<"Errore in CategoricalCluster.ViewThetaRightId"<<std::endl;
        std::cerr<<"Id eccede rispetto gli elemtni contenuti"<<std::endl;
        exit(1);
    }
    return ThetaRight[_id];
}

//Estre il peso globale del cluster
double CategoricalCluster::ViewBeta() const{

    return Beta;
}

//Estrae il peso globale del sub-cluster sinistro
double CategoricalCluster::ViewBetaLeft() const{

    return BetaLeft;
}

//Estrae il peso globale del sub-cluster  destro
double CategoricalCluster::ViewBetaRight() const{

    return BetaRight;
}

//Estrae i valori delle statistiche del cluster
void CategoricalCluster::ViewStatistics( STAT& _c) const{    /////////

    _c=c;
}

//Estrae i valori delle statistiche del sub-cluster sinistro
void CategoricalCluster::ViewStatisticsLeft( STAT& _cLeft) const{ //////////

    _cLeft=cLeft;
}

//Estrae i valori delle statistiche del sub-cluster destro
void CategoricalCluster::ViewStatisticsRight( STAT& _cRight) const{ /////////

    _cRight=cRight;
}

//Estrae il numero di tavoli globali nel cluster
unsigned int CategoricalCluster::ViewGlobalTable()const{

    return NrTable;
}

//Estrae il numero di tavoli globali nel sub-cluster sinistro
unsigned int CategoricalCluster::ViewGlobalTableLeft() const{

    return NrTableLeft;
}

//Estrae il numero di tavoli globali nel sub-cluster destro
unsigned int CategoricalCluster::ViewGlobalTableRight() const{

    return NrTableRight;
}

// Azzera le statistiche nel cluster
void CategoricalCluster::ResetStatistics(unsigned int W){

	c.assign(W,0);
}

//Azzera le statiche nel sub-cluster sinistro
void CategoricalCluster::ResetStatisticsLeft(unsigned int W){

	cLeft.assign(W,0);
}

//Azzera le statistiche sel sub-cluster destro
void CategoricalCluster::ResetStatisticsRight(unsigned int W){

	cRight.assign(W,0);
}

//Aggiorna le statistiche, ovvero gli iperparametri dei parametri latenti dei cluster e sub-cluster
void CategoricalCluster::UpdateStatistics(STAT& counts4cleft, STAT& counts4cright ){

	if(counts4cleft.size() != counts4cright.size()){
		std::cerr<<" error in UpdateStatistics: counts4cleft and counts4cleft must have the same dimension "<<std::endl;
		exit(1);
	}

	if(cLeft.empty())
		cLeft.assign(counts4cleft.size(),0);

	if(cRight.empty())
		cRight.assign(counts4cright.size(),0);

	if(c.empty())
		c.assign(cLeft.size(),0);

	for(typename STAT::size_type w = 0; w<counts4cleft.size(); ++w){
		cLeft[w] += counts4cleft[w];
		cRight[w] += counts4cright[w];
		c[w] += counts4cleft[w] + counts4cright[w];

	}

}

//Controlla se il cluster non contiene elementi
bool CategoricalCluster::IsEmpty() const {

    unsigned int sum = 0;
    for(auto i: c)
        sum += i;

    if(sum == 0)
        return true;
    else
        return false;
}

#endif
