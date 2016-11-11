#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <vector>
#include <iostream>
#include "Type.hpp"

using std::vector;

/*! \file Cluster.hpp
 *
 * \brief Strutture dei dati per la gestione dei cluster in base al modello scelto.
 * In queste classi si definisce solo come gestire i parametri lanteti e gli altri parametri del cluster.
 * Non sono presenti campionamenti, ma solo metodi che stampano e fissano i valori dei cluster e sub-cluster
 * \date Febbraio 2016
 */


/*! \brief Modello Generico di Cluster
 *
 *  Classe astratta dove tutti i metodi virtuali sono null.
 *  Le classi che ereditano da Cluster Generic servono per estrarre, memorizzare ed impostare
 *  i dati associati ai cluster ed ai sub-cluster. Pertanto si tratta di classi di appoggio, non fanno nessun tipo di campionamento
 * \date Febbraio 2016
 */

template <template <unsigned int> class ClassType, unsigned int DIM>
class GenericCluster
{
public:


    /*! \brief  Fissa il parametro latente del cluster
     *  \param   _Theta - parametro latente in ingresso di tipo THETA
     */
    virtual void SetTheta(typename ClassType<DIM>::THETA&) = 0;

    /*! \brief  Fissa il parametro latente del sub-cluster sinistro
     *  \param  _ThetaLeft - parametro latente in ingresso di tipo THETA
     */
    virtual void SetThetaLeft(typename ClassType<DIM>::THETA&) = 0;

    /*! \brief  Fissa il parametro latente del sub-cluster destro
     *  \param  _ThetaRight - parametro latente in ingresso di tipo THETA
     */
    virtual void SetThetaRight(typename ClassType<DIM>::THETA&) = 0;

    /*! \brief  Fissa il peso del cluster
     *  \param  _Beta - peso in ingresso
     */
    virtual void SetBeta(double) = 0;

    /*! \brief	Fissa il peso del sub-cluster sinistro
     *  \param _BetaLeft - peso in ingresso
     */
    virtual void SetBetaLeft(double) = 0;

    /*! \brief	Fissa il peso del sub-cluster destro
     *  \param _BetaRight - peso in ingresso
     */
    virtual void SetBetaRight(double) = 0;

    /*! \brief Fissa il numero di tavoli globali che caratterizza il cluster
        \param NrTable - numero di tavoli che caratterizza il cluster
	 */
    virtual void SetGlobalTable (unsigned int) = 0;

    /*! \brief Fissa il numero di tavoli globali che caratterizza il sub-cluster sinistro
        \param NrTableLeft - numero di tavoli che caratterizza il sub-cluster sinistro
	 */
    virtual void SetGlobalTableLeft (unsigned int) = 0;

    /*! \brief Fissa il numero di tavoli globali che caratterizza il sub-cluster destro
        \param NrTableRight - numero di tavoli che caratterizza il sub-cluster destro
	 */
    virtual void SetGlobalTableRight (unsigned int) = 0;

	/*! \brief  Fissa le statistiche, ovvero gli iperparametri dei parametri latenti del cluster
     *  \param  _c - statistiche del cluster
	 */
	virtual void SetStatistics(typename ClassType<DIM>::STAT&) = 0;

    /*! \brief  Fissa le statistiche, ovvero gli iperparametri dei parametri latenti del sub-cluster sinistro
     *  \param  _cLeft - statistiche del sub-cluster sinistro
	 */
    virtual void SetStatisticsLeft(typename ClassType<DIM>::STAT&) = 0;

    /*! \brief  Fissa le statistiche, ovvero gli iperparametri dei parametri latenti del sub-cluster destro
     *  \param _cRight - statistiche del sub-cluster destro
	 */
    virtual void SetStatisticsRight(typename ClassType<DIM>::STAT&) = 0;

    /*! \brief Estrae il parametro latente del cluster
        \param _Theta - oggetto di tipo THETA in cui viene memorizzato il parametro latente del cluster
     */
    virtual void ViewTheta (typename ClassType<DIM>::THETA&) const = 0;

    /*! \brief Estrae il parametro latente del sub-cluster sinistro
        \param _ThetaLeft - oggetto di tipo THETA in cui viene memorizzato il parametro latente del sub-cluster sinistro
     */
    virtual void ViewThetaLeft(typename ClassType<DIM>::THETA&) const = 0;

    /*! \brief Estrae il parametro latente del sub-cluster destro
        \param _ThetaRight - oggetto di tipo THETA in cui viene memorizzato il parametro latente del sub-cluster destro
     */
    virtual void ViewThetaRight(typename ClassType<DIM>::THETA&) const = 0;

    /*! \brief Estra il peso globale del cluster
     *  \return Peso del cluster
     */
    virtual double ViewBeta() const = 0;

    /*! \brief Estra il peso globale del sub-cluster sinistro
     *  \return Peso del sub-cluster sinistro
     */
    virtual double ViewBetaLeft() const = 0;

    /*! \brief Estra il peso globale del sub-cluster destro
     *  \return Peso del del sub-cluster destro
     */
    virtual double ViewBetaRight() const = 0;

    /*! \brief Estra la statistica del parametro latente del cluster
     *  \param _c - oggetto di tipo STAT in cui viene memorizzato la statistica del parametro latente del cluster
     */
    virtual void ViewStatistics( typename ClassType<DIM>::STAT&) const = 0;

    /*! \brief Estra la statistica del parametro latente del sub-cluster sinistro
     *  \param _cLeft - oggetto di tipo STAT in cui viene memorizzato la statistica del parametro latente del sub-cluster sinistro
     */
    virtual void ViewStatisticsLeft( typename ClassType<DIM>::STAT&) const = 0;

    /*! \brief Estra la statistica del parametro latente del sub-cluster destro
     *  \param _cRight - oggetto di tipo STAT in cui viene memorizzato la statistica del parametro latente del sub-cluster destro
     */
    virtual void ViewStatisticsRight( typename ClassType<DIM>::STAT&) const = 0;

    /*! \brief Estra il numero globale dei tavoli nel cluster
     *  \return Numero di tavoli nel cluster
     */
    virtual unsigned int ViewGlobalTable () const = 0;

    /*! \brief Estra il numero globale dei tavoli nel sub-cluster sinistro
     *  \return Numero di tavoli nel sub-cluster sinistro
     */
    virtual unsigned int ViewGlobalTableLeft () const = 0;

    /*! \brief Estra il numero globale dei tavoli nel sub-cluster destro
     *  \return Numero di tavoli nel sub-cluster destr
     */
    virtual unsigned int ViewGlobalTableRight () const = 0;

    /*! \brief  Azzera le statistiche nel cluster
     *  \param  W - dimensione della statistica con cui aggirnare gli iperparametri del parametro latente del cluster
	 */
	virtual void ResetStatistics(unsigned int) = 0;

    /*! \brief  Azzera le statistiche nel sub-cluster sinistro
     *  \param  W - dimensione della statistica con cui aggirnare gli iperparametri del parametro latente del sub-cluster sinistro
	 */
	virtual void ResetStatisticsLeft(unsigned int) = 0;

	/*! \brief  Azzera le statistiche nel sub-cluster destro
     *  \param  W - dimensione della statistica con cui aggirnare gli iperparametri del parametro latente del sub-cluster destro
	 */
	virtual void ResetStatisticsRight(unsigned int) = 0;

	/*! \brief  Aggiorna le statistiche, ovvero gli iperparametri dei parametri latenti dei cluster e sub-cluster
     *  \param  counts4cleft - statistiche per aggiornare gli iperparametri dei parametri latenti del sub-cluster sinistro
     *  \param  counts4cright - statistiche per aggiornare gli iperparametri dei parametri latenti del sub-cluster destro
	 */
	virtual void UpdateStatistics(typename ClassType<DIM>::STAT&, typename ClassType<DIM>::STAT&) = 0;

    /*! \brief Controlla se il cluster non contiene elementi
     *  \return TRUE se il cluster e' vuoto FALSE se non lo e'
     */
    virtual bool IsEmpty() const = 0;

};

/*! \brief Gestione informazioni cluster e sub-cluster per dati di verosimiglianza Categorical
 *
 *  Questa classe si occupa di memorizzare ed estrarre informazioni riguardati il peso globale del cluster e sub-cluster, parametri latenti
 *  dei cluster e sub-cluster, informazioni inerenti agli aggiornamenti degli iperparametri dei parametri lantenti. Nel caso di
 *  verosimiglianza caregorica. i dati sono ripetuti, i parametri lantenti sono i pesi degli elementi distinti e le statistiche
 *  sono i conteggi degli elementi distinti nel cluster.
 *  I parametri latenti non sono altro che i parametri della mistura.
 *  \date Febbrario 2016
 */
class CategoricalCluster final: GenericCluster<TypeCategorical,1>
{
    public:

        //DEFINIZIONE DEI TIPI PUBBLICI

        /*! \brief Parametro latente, vettore di pesi degli elementi distinti nel cluster
         */
        using THETA = TypeCategorical<1>::THETA;

        /*! \brief Singlo dato, dati ripetuti
         */
        using Point = TypeCategorical<1>::Point;

        /*! \brief Statistiche per aggiornare l'iperparametro del parametro latente, numero di dati che sono cotentuti nel cluster e sub-cluster
         */
        using STAT  = TypeCategorical<1>::STAT;

    private:

        // DEFINIZIONE TIPI E STRUTTURE DATI PRIVATI

        /*! \brief Peso globale del cluster
         */
        double Beta;

        /*! \brief Peso globale del sub-cluster sinistro
         */
        double BetaLeft;

         /*! \brief Peso globale del sub-cluster destro
         */
        double BetaRight;

        /*! \brief Parametro latente del cluster: peso degli elementi distinti nel cluster
         */
        THETA Theta;

        /*! \brief Parametro latente del sub-cluster sinisto: peso degli elementi distinti nel sub-cluster sinistro
         */
        THETA ThetaLeft;

        /*! \brief Parametro latente del sub-cluster destro: peso degli elementi distinti nel sub-cluster destro
         */
        THETA ThetaRight;

        /*! \brief Statisitca per aggiornare il parametro latente del cluster: conteggi degli elementi finiti nel cluster
         */
        STAT c;

        /*! \brief Statisitca per aggiornare il parametro latente del sub-cluster sinistro: conteggi degli elementi finiti nel sub-cluster sinistro
         */
        STAT cLeft;

        /*! \brief Statisitca per aggiornare il parametro latente del sub-cluster destro: conteggi degli elementi finiti nel sub-cluster destro
         */
        STAT cRight;

        /*! \brief Numero di tavoli nel cluster
         */
        unsigned int NrTable;

         /*! \brief Numero di tavoli nel sub-cluster sinistro
         */
        unsigned int NrTableLeft;

         /*! \brief Numero di tavoli nel sub-cluster destro
         */
        unsigned int NrTableRight;

    public:

        //COSTRUTTURI E DISTRUTTORE

        /*! \brief Costruttore di default
         */
        CategoricalCluster();

        /*! \brief Distruttore di default
         */
        ~CategoricalCluster()=default;

        /*! \brief Costruttore che richiede tutte le informazioni del cluster e subcluster
         *  \param _Beta - Peso globale del cluster
         *  \param _BetaLeft - Peso globale del sub-cluster sinistro
         *  \param _BetaRight - Peso globale del sub-cluster destro
         *  \param _Theta - Parametro latente del cluster: peso degli elementi distinti nel cluster
         *  \param _ThetaLeft - Parametro latente del sub-cluster sinistro: peso degli elementi distinti nel sub-cluster sinistro
         *  \param _ThetaRoght - Parametro latente del sub-cluster destro: peso degli elementi distinti nel sub-cluster destro
         *  \param _c - Statisitca per aggiornare il parametro latente del cluster: conteggi degli elementi finiti nel cluster
         *  \param _cLeft - Statisitca per aggiornare il parametro latente del sub-cluster sinistro: conteggi degli elementi finiti nel sub-cluster sinistro
         *  \param _cRight - Statisitca per aggiornare il parametro latente del sub-cluster destro: conteggi degli elementi finiti nel sub-cluster destro
         *  \param -NrTable - Numero di tavoli nel cluster
         *  \param -NrTableLeft - Numero di tavoli nel sub-cluster sinistro
         *  \param -NrTableRight - Numero di tavoli nel sub-cluster destro
         */
        CategoricalCluster(double _Beta, double _BetaLeft, double _BetaRight, THETA& _Theta ,THETA& _ThetaLeft, THETA& _ThetaRight, STAT& _c, STAT& _cLeft, STAT& _cRight, unsigned int _NrTable, unsigned int _NrTableLeft, unsigned int _NrTableRight);

       //METODI PER FISSARE
        /*! \brief  Fissa il parametro latente del cluster
         *  \param   _Theta - parametro latente in ingresso di tipo THETA
         */
        void SetTheta(THETA& _Theta);

        /*! \brief  Fissa il parametro latente del sub-cluster sinistro
         *  \param  _ThetaLeft - parametro latente in ingresso di tipo THETA
         */
        void SetThetaLeft(THETA& _ThetaLeft);

        /*! \brief  Fissa il parametro latente del sub-cluster destro
         *  \param  _ThetaRight - parametro latente in ingresso di tipo THETA
         */
        void SetThetaRight(THETA& _ThetaRight);

        /*! \brief  Fissa il peso del cluster
         *  \param  _Beta - peso in ingresso
         */
        void SetBeta(double _Beta);

        /*! \brief	Fissa il peso del sub-cluster sinistro
         *  \param _BetaLeft - peso in ingresso
         */
        void SetBetaLeft(double _BetaLeft);

        /*! \brief	Fissa il peso del sub-cluster destro
         *  \param _BetaLeft - peso in ingresso
         */
        void SetBetaRight(double _BetaRight);

        /*! \brief Fissa il numero di tavoli globali che caratterizza il cluster
         *  \param NrTable - numero di tavoli che caratterizza il cluster
         */
        void SetGlobalTable(unsigned int _NrTable);

        /*! \brief Fissa il numero di tavoli globali che caratterizza il sub-cluster sinistro
         *  \param NrTableLeft - numero di tavoli che caratterizza il sub-cluster sinistro
         */
        void SetGlobalTableLeft(unsigned int _NrTableLeft);

        /*! \brief Fissa il numero di tavoli globali che caratterizza il sub-cluster destro
         *  \param NrTableLeft - numero di tavoli che caratterizza il sub-cluster destro
         */
        void SetGlobalTableRight(unsigned int _NrTableRight);

        /*! \brief  Fissa le statistiche, ovvero gli iperparametri dei parametri latenti del cluster
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

        /*! \brief Estra la statistica del parametro latente del cluster
         *  \param _c - oggetto di tipo STAT in cui viene memorizzato la statistica del parametro latente del cluster
         */
        void ViewStatistics( STAT& _c) const;

        /*! \brief Estra la statistica del parametro latente del sub-cluster sinistro
         *  \param _cLeft - oggetto di tipo STAT in cui viene memorizzato la statistica del parametro latente del sub-cluster sinistro
         */
        void ViewStatisticsLeft( STAT& _cLeft) const;

        /*! \brief Estra la statistica del parametro latente del sub-cluster destro
         *  \param _cLeft - oggetto di tipo STAT in cui viene memorizzato la statistica del parametro latente del sub-cluster destro
         */
        void ViewStatisticsRight( STAT& _cRight) const;

        /*! \brief  Azzera le statistiche nel cluster
         *  \param  W - dimensione della statistica con cui aggirnare gli iperparametri del parametro latente del cluster, numero di elementi distinti
         */
        void ResetStatistics(unsigned int W);

        /*! \brief  Azzera le statistiche nel sub-cluster sinistro
         *  \param  W - dimensione della statistica con cui aggirnare gli iperparametri del parametro latente del sub-cluster sinistro,
                    numero di elementi distinti
         */
        void ResetStatisticsLeft(unsigned int W);

        /*! \brief  Azzera le statistiche nel sub-cluster destro
         *  \param  W - dimensione della statistica con cui aggirnare gli iperparametri del parametro latente del sub-cluster destro,
                    numero di elementi distinti
         */
        void ResetStatisticsRight(unsigned int W);

        /*! \brief  Aggiorna le statistiche, ovvero gli iperparametri dei parametri latenti dei cluster e sub-cluster
         *  \param  counts4cleft - statistiche per aggiornare gli iperparametri dei parametri latenti del sub-cluster sinistro
         *  \param  counts4cright - statistiche per aggiornare gli iperparametri dei parametri latenti del sub-cluster destro
         */
        void UpdateStatistics(STAT& counts4cleft, STAT& counts4right);

        /*! \brief Controlla se il cluster non contiene elementi
         *  \return TRUE se il cluster e' vuoto FALSE se non lo e'
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
