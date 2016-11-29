#ifndef __STRUCT__HPP__
#define __STRUCT__HPP__

#include <vector>
#include <utility>
#include <unordered_map>
#include "Type.hpp"

using std::vector;
using std::pair;
using std::unordered_map;

/*! \file Struct.hpp
 * \brief Gathers structures used in the methods of HDP_MCMC class.
 * We chose to create a separate file for the structures' definition because they are common to different methods.
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */

using POINT = TypeCategorical<1>::Point;
using STAT  = TypeCategorical<1>::STAT;


/*! \brief Clusters' global weights
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */
struct BETA {

    double a;
    vector<double> b_c;
	//servono nel Gibbs_SubTopic
    vector<double> Left;
    vector<double> Right;

	double k;

};

/*! \brief Number of elements of group j in cluster k
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */
struct NJK {

    unsigned int a;
    unsigned int b;
    unsigned int c;

	unsigned int k;

    pair<unsigned int, unsigned int> a_sub;
    pair<unsigned int, unsigned int> b_sub;
    pair<unsigned int, unsigned int> c_sub;

};

/*! \brief Group specific clusters' weights
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */
struct PI {

    vector<double> b_c;
    double a;

    vector<double> Tilde_b_c;
    //proposte dopo il campionamento dei subTopic
    double Left;
    double Right;

};

/*! \brief Tables
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */
struct NUMTABLE {

    unsigned int a_Left;
    unsigned int a_Right;
    unsigned int a;
	unsigned int ja_Left;
    unsigned int ja_Right;
    unsigned int ja;
    unsigned int Tilde_sum;

    //sarebbero gli m_tildaLeft e m_tildaRight
    vector<unsigned int> Tilde_b_c;
	vector<unsigned int> Tilde_k;
};

/*! \brief Statistics
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */
struct C {
	STAT a; //conterrà c_a_id per ogni id da 1 a W
	STAT b;   // estraggo i conteggi dai cluster esistenti
	STAT c;
	STAT a_left; //conterrà c_a_id per ogni id da 1 a W
	STAT b_left;   // estraggo i conteggi dai cluster esistenti
	STAT c_left;
	STAT a_right; //conterrà c_a_id per ogni id da 1 a W
	STAT b_right;   // estraggo i conteggi dai cluster esistenti
    STAT c_right;

};

/*! \brief Structure for a cluster
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */
struct CLUSTER{

	pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>> b;
	pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>> c;
	pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>> a; //da costruire con le parti dx e sx
	unordered_map<POINT,unsigned int> a_sx; //parte sx del cluster a, da riempire per ogni documento
	unordered_map<POINT,unsigned int> a_dx;
};

/*! \brief  Structure for data counts
 * \authors{Debora Parisi and Stefania Perego}
 * \date February 2016
 */
struct DATACOUNT{

	vector<unsigned int> b;   //conterrà njb per ogni j
	vector<unsigned int> c;   //conterrà njc per ogni j
	vector<unsigned int> a;   //conterrà nja_cap per ogni j

};

#endif
