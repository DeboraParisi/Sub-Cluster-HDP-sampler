
#include "../src/include/GetPot"
#include "CategoricalDatasetGenerator.cpp"
#include <iostream>
//void CategoricalDatasetGenerator(unsigned int K, double Gamma, unsigned int D, double Alpha , unsigned int W, double Lambda, const vector<unsigned int>& Nj);

int main(int argc, char** argv){

    GetPot command_line (argc,argv);
    
    std::string DataFile;
    
    if(!command_line.search("-f"))
    	DataFile = command_line.follow(" ",1," ");
    else
	DataFile = command_line.follow("GeneratorDataFile.txt",1,"-f");
	
	std::cout<<DataFile<<std::endl;
    GetPot file (DataFile.c_str());
	
	//bool is_K = static_cast<bool>(file("is_K", 1));
	unsigned int K = file("K", 5);
	
	//bool is_gamma = static_cast<bool>(file("is_gamma", 1));
	double gamma = file("gamma",10.0);
	
	//bool is_D = static_cast<bool>(file("is_D", 1));
	unsigned int D = file("D", 5);
	
	//bool is_alpha = static_cast<bool>(file("is_alpha", 1));
	double alpha = file("alpha",1.0);
	
	//bool is_W = static_cast<bool>(file("is_W", 1));
	unsigned int W = file("W", 100);
	
	//bool is_lambda = static_cast<bool>(file("is_lambda", 1));
	double lambda = file("lambda",0.5);
	
	//bool is_Nj = static_cast<bool>(file("is_Nj", 1));
	unsigned int Nj = file("Nj",50);

    vector<unsigned int> AllNj(D,Nj);

    CategoricalDatasetGenerator(K,gamma,D,alpha,W,lambda,AllNj);

    return 0;
}
