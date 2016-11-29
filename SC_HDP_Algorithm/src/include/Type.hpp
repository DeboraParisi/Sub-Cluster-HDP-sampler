#ifndef _TYPE_HPP_
#define _TYPE_HPP_

#include <vector>
#include <iostream>

using std::vector;


/*! \file Type.hpp
 *
 * \brief  Data's structures which define the model's type and clusters' type
 */

/*! \brief This class defines data's types, with them it is possible to represent the categorcal likelihood
 *
 *	\authors { Debora Parisi and Perego Stefania }
 *  \date February 2016
*/
template < unsigned int DIM=1>
class TypeCategorical
{
public:

    /*! \brief Vector of latent parameter of cluster or sub-cluster
    */
    using THETA = vector<double> ;
    /*!  \brief Type of single observation: Xji
    */
    using Point = unsigned int ;
    /*! \brief Vector of latent parameter's hyperparameter of cluster or sub-cluster
    */
	using HYP = vector<double>;
	/*! \brief Statistics' vector. In this container there are the update of latent parameters' hyperparameters of clusters or sub-clusters.
	 */
	using STAT = vector<unsigned int>;

};

#endif
