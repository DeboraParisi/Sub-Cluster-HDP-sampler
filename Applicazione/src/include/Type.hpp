#ifndef _TYPE_HPP_
#define _TYPE_HPP_

#include <vector>
#include <iostream>

using std::vector;

/*! \file Type.hpp
 *
 * \brief  Strutture dati per i tipi di modello e dei cluster
 */

/*! \brief Classe dei tipi per dati con verosimiglianza categorica
 *
 *  \date Febbraio 2016
*/
template < unsigned int DIM=1>
class TypeCategorical
{
public:

    /*! \brief Vettore dei parametri latenti dei cluster o sub-cluster
    */
    using THETA = vector<double> ;
    /*!  \brief Dato singolo
    */
    using Point = unsigned int ;
    /*! \brief Vettore degli iperparametri dei parametri lantenti dei cluster o sub-cluster
    */
	using HYP = vector<double>;
	/*! \brief Vettore delle statistiche, aggiornamenti degli iperparametri dei parametri latenti dei cluster o sub-cluster
	*/
	using STAT = vector<unsigned int>; //conteggi per aggiornare iperparametri

};

#endif
