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
 *
 * \brief Raccolta delle strutture che vengono utilizzate nei metodi di HDP_MCMC.hpp
 * Si e' scelto di creare un file a parte con la dichiarazione delle strutture, perchè il loro
 * utilizzo e' comune a più metodi
 *
 * Ogni volta che si aggiunge un nuovo modello e dei nuovi tipi di dati, bisogna aggiornare la lista in cui si rinomina i tipi di dati
 *
 * \date Febbraio 2016
 */

using POINT = TypeCategorical<1>::Point;
using STAT  = TypeCategorical<1>::STAT;


/*! \brief Pesi globali dei cluster
 */
struct BETA {

    double a;
    vector<double> b_c;
	//servono nel Gibbs_SubTopic
    vector<double> Left;
    vector<double> Right;

	double k;

};

/*! \brief Numero degli elementi del gruppo j che sono nel cluster k
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

/*! \brief Pesi dei cluster in ogni gruppo
 */
struct PI {

    vector<double> b_c;
    double a;

    vector<double> Tilde_b_c;
    //proposte dopo il campionamento dei subTopic
    double Left;
    double Right;

};

/*! \brief Tavoli
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

/*! \brief Statistiche
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

struct CLUSTER{

	pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>> b;
	pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>> c;
	pair<unordered_map<POINT,unsigned int>,unordered_map<POINT,unsigned int>> a; //da costruire con le parti dx e sx
	unordered_map<POINT,unsigned int> a_sx; //parte sx del cluster a, da riempire per ogni documento
	unordered_map<POINT,unsigned int> a_dx;
};

struct DATACOUNT{

	vector<unsigned int> b;   //conterrà njb per ogni j
	vector<unsigned int> c;   //conterrà njc per ogni j
	vector<unsigned int> a;   //conterrà nja_cap per ogni j

};

#endif
