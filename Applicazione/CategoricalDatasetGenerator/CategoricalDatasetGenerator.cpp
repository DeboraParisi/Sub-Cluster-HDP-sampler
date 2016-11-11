
#include "../src/include/omprng.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

using std::vector;
using std::unordered_map;


void CategoricalDatasetGenerator(unsigned int K, double Gamma, unsigned int D, double Alpha , unsigned int W, double Lambda, const vector<unsigned int>& Nj);


void CategoricalDatasetGenerator(unsigned int K, double Gamma, unsigned int D, double Alpha , unsigned int W, double Lambda, const vector<unsigned int>& Nj){

   	std::ofstream file("Dataset.txt");
	
	std::ofstream file2("Variables.txt");

	std::ofstream file3 ("Vocabulary.txt");
	
	if (file.fail()){
		std::cerr <<"ERROR GENERATED (CategoricalDatasetGenerator)!" <<std::endl;
		std::cerr <<"Cannot open Dataset.txt" <<std::endl;
		exit(1);
	}
	
	if (file2.fail()){
		std::cerr <<"ERROR GENERATED (CategoricalDatasetGenerator)!" <<std::endl;
		std::cerr <<"Cannot open Variables.txt" <<std::endl;
		exit(1);
	}

	if (file3.fail()){
		std::cerr <<"ERROR GENERATED (CategoricalDatasetGenerator)!" <<std::endl;
		std::cerr <<"Cannot open Vocabulary.txt" <<std::endl;
		exit(1);
	}

    if( Nj.size() != D ){
	  std::cerr<< " In CategoricalDatasetGenerator, Nj must have dimension D"<<std::endl;
	  exit(1);
	}

    omprng Gen;
	Gen.setNumThreads(1);

    vector<double> Beta_Sample;
    Beta_Sample.reserve(K);
	
	for(unsigned int k=0; k<K; ++k)
	    Beta_Sample.push_back(Gen.rbeta(1,Gamma));  
	
	vector<double> Beta_Sample_tilde;
    Beta_Sample_tilde.reserve(K);
	
	for(unsigned int k=0; k<K; ++k)
	    Beta_Sample_tilde.push_back(1-Beta_Sample[k]); 
		
	// calcolo i prodotti cumulati

    vector<double> cumprod;	
	cumprod.resize(K);
	
	cumprod[0]=1;
	
	double prod = 1;
	
	for(unsigned int k=1; k<K; ++k){
	   for(unsigned int l=0; l<k; ++l)
	       prod *= Beta_Sample_tilde[l];
	   cumprod[k] = prod;
	   prod = 1;
	   
	}   
	
	
	// definisco il vettore dei pesi definitivo
    vector<double> Beta;
    Beta.resize(K);
	//double betasum = 0.0;
	//double checksum = 0.0;
	
	for(unsigned int k=0; k<K; ++k){
	   Beta[k] = Beta_Sample[k]*cumprod[k];
	   //betasum += Beta[k];
	}

	
	unordered_map<unsigned int,vector<double>> PiContainer;
	vector<double> PIj_sample;
	vector<double> PIj_sample_tilde;
	vector<double> PIj;
	double sum = 0.0;
	//double pisum = 0.0;
	
	for(unsigned int d=0; d<D; ++d){
	   
	   PIj_sample.reserve(K);
	   PIj.resize(K);
	   for(unsigned int k=0; k<K; ++k){
	       for(unsigned int l=0; l<k+1; l++)
	           sum+= Beta[l];
	       PIj_sample.push_back ( Gen.rbeta(Alpha*Beta[k], Alpha*(1-sum )));
		   //std::cout<< "Alpha*Beta[k] "<< Alpha*Beta[k]<<" Alpha*(1-sum ) "<<Alpha*(1-sum )<<std::endl;
	       sum = 0.0;
	   } 
	   
	
    PIj_sample_tilde.reserve(K);
	
	for(auto i: PIj_sample)
	   PIj_sample_tilde.push_back(1-i);

    	// calcolo i prodotti cumulati	
	cumprod.clear();	
	cumprod.resize(K);
	cumprod[0]=1;
	prod = 1;
	
	for(unsigned int k=1; k<K; ++k){
	   for(unsigned int l=0; l<k; ++l)
	       prod *= PIj_sample_tilde[l];
	   cumprod[k] = prod;
	   prod = 1;	   
	} 
	 
	//vettore PIj definitivo
	
	//sum = 0.0;
    	
    for(unsigned int k=0; k<K; ++k){
	   PIj[k] = PIj_sample[k]*cumprod[k];	 
	  // pisum += PIj[k];
	} 

    PiContainer.insert({d,PIj});	
	

    PIj_sample.clear();
    PIj.clear();
	PIj_sample_tilde.clear();
	//pisum = 0.0;
    }	
	
	
	unordered_map< unsigned int, vector<double> > Theta;
    vector<double> Params_Theta (W,Lambda);
    vector<double> Temp_Theta;

    for(unsigned int k=0; k<K; k++){
	
        Temp_Theta.clear();
        Temp_Theta.resize(W);

        Gen.rdirichlet(Params_Theta,Temp_Theta);

        Theta.insert({k,Temp_Theta});
		
    }
	
	// Gen.rdiscrete normalizza i pesi, se non lo sono già
	
	unsigned int z;
	unsigned int id;
	unordered_map< unsigned int, unsigned int > id_count;  // ogni documento ha la sua
	
	
	for(unsigned int d=0; d<D; ++d){
	    
		for(unsigned int n=0; n< Nj[d] ; ++n){
	
	       z = Gen.rdiscrete( PiContainer[d]);  
		   id = Gen.rdiscrete(Theta[z]);
		//std::cout<< "d "<<d<<std::endl;
	    //std::cout<< "z "<<z<<std::endl;
        //std::cout<< "id "<<id<<std::endl;	
		
           if(!(id_count.insert({id,1}).second))
              ++ id_count[id]; 	   
	    
		}
		
		for(auto iter : id_count)
		   file << d+1 << " " << (iter.first + 1) << " " << iter.second <<std::endl;
		   
		id_count.clear();   
	
	}
	
	// ora credo il file con le variabili D,W,N
	
	unsigned int N = 0;
	
	for(auto i: Nj)
	   N += i;
	
	file2<<D<<std::endl;
	file2<<W<<std::endl;
	file2<<N<<std::endl;		   

   	file.close();
	file2.close();

	//creo il vocabolario
	for(unsigned int i=1; i<=W; i++)
		file3<<i<<endl;
	file3.close();

}







