#include "../src/include/PosteriorAnalysis.hpp"
#include "../src/include/GetPot"
#include "../src/include/Type.hpp"
#include <RInside.h>
#include <iostream>
#include <vector>
#include <string>


int main(int argc, char* argv[]){

	using std::cout;
	using std::cin;
	using std::endl;
	using std::cerr;

	GetPot command_line (argc,argv);
    const std::string DataFile = command_line.follow("DataFile.txt",1,"-f");

    GetPot file (DataFile.c_str());

    const std::string wd="setwd(\"~/APSC_Parisi_Perego/PosteriorAnalysis/R_results\")";		// imposta la wd per R, la passo ad ogni funzione che usa R

	unsigned int W;   //  nr parole distinte
    unsigned int D;   // nr documenti
	unsigned int N;   // numerosità dataset

	bool is_LPML = static_cast<bool> (file ("is_LPML", 0));
	const unsigned long burnin = file("burnin",0);
	bool is_alpha_prior = static_cast<bool>(file("is_alpha_prior",0));
    bool is_gamma_prior = static_cast<bool>(file("is_gamma_prior",0));
	const unsigned long Iterations = file("Iter",10000);
    unsigned int Kdim = 0;
	unsigned int Alphadim = 0;
	unsigned int Gammadim = 0;
	char AlphaTry = 'y';
	char GammaTry = 'y';
	unsigned int AlphaBurnin = 0;
	unsigned int GammaBurnin = 0;
	unsigned int AlphaThinning = 1;
	unsigned int GammaThinning = 1;
	char Try = 'y';
	unsigned int Burnin = 0;
	unsigned int Thinning = 1;
	char Prior;

	std::ifstream var_file("Variables.txt");
    vector<unsigned int> Variable;

    if(var_file.fail()){
        std::cerr<<"ERROR in Variables acquisition"<<std::endl;
        std::cerr<<"I cannot open Variables.txt"<<std::endl;
        exit(1);
    }

    std::string ss1;
    unsigned int temp;

    while(std::getline(var_file,ss1)){
        std::istringstream SSTR(ss1);
        SSTR>>temp;
        Variable.push_back(temp);
    }

    if(Variable.size()!=3){
        std::cerr<<"ERROR in Variables acquisition"<<std::endl;
        std::cerr<<"Less data than I expect"<<std::endl;
        exit(1);
    }

    D=Variable[0];
    W=Variable[1];
    N=Variable[2];

	cout<<"Starting Posterior Analysis..."<<endl;
    cout<<endl;

    // Occorre un oggetto RInside che viene passato per referenza alla funzione
    cout <<"Opening R environment... "<<endl;
    RInside R(argc, argv);
    cout <<"Done. " <<endl;
	cout<<endl;

    cout<<"Istantiating CategoricalPosteriorAnalysis object..."<<endl;
	CategoricalPosteriorAnalysis<1> CPA;
	CPA.SetW(W);
	CPA.SetD(D);
	CPA.SetN(N);
	cout<<"Done"<<endl;
	cout<<endl;

	cout<<"Setting working directory..."<<endl;
	CPA.Setwd(wd);
	cout<<"Done"<<endl;
	cout<<endl;
	
	cout<<"Loading all K..."<<endl;
	Kdim = CPA.SetAllK();
	cout<<"Done."<<endl;
	cout<<endl;
	
	cout<< "Loading theta and beta..."<<endl;
    CPA.SetTheta(Iterations);
	cout<<"Done"<<endl;
	cout<<endl;
	
	cout<< "Loading pi..."<<endl;
    CPA.SetPi(Iterations);
	cout<<"Done"<<endl;
	cout<<endl;
	
	cout<< "Loading docs names..."<<endl;
	CPA.SetDocs();
	cout<<"Done"<<endl;
	cout<<endl;
	
	
	unsigned long BestClustering = 0;

	//if(N<10001){
		cout<<"Loading labels..."<<endl;
		//if(Iterations >= 10000)
			//CPA.LoadLabels(10000);
		//else
			CPA.LoadLabels(Iterations);
		cout<<"Done"<<endl;
		cout<<endl;

		cout<<"Finding least square clustering..."<<endl;
		cout<<endl;
		BestClustering = CPA.LeastSquareClustering();
		cout<<"Done"<<endl;
		cout<<endl;
    //}

		cout<<"Associating docs..."<<endl;
		CPA.AssociatingDocs();
		cout<<"Done."<<endl;
		cout<<endl;
	
	
	/*
	if(is_alpha_prior && is_gamma_prior){

		cout<<"Would you like to do the Posterior Analysis on Alpha and Gamma? Possible only if they aren't fixed"<<endl;
		cout<<"Press y [Yes] or n [no]  ";
		cin>>Prior;
		while(Prior !='y' && Prior != 'n'){
			cout<<"press y [yes] or n [no] "<<endl;
			cin>>Prior;
		}

		cout<<endl;
	}

	if(Prior == 'y'){
		cout<<"Loading Alpha..."<<endl;
		Alphadim = CPA.SetAllAlpha();
		cout<<"Done."<<endl;
		cout<<endl;

		cout<<"Loading Gamma..."<<endl;
		Gammadim = CPA.SetAllGamma();
		cout<<"Done."<<endl;
		cout<<endl;
	}
	*/

	/*cout<<"Loading vocabulary..."<<endl;
	CPA.SetVocabulary();
	cout<<"Done"<<endl;
	cout<<endl;

	cout<<"Tracking Clusters in the last 100 iterations of the HDP_MCMC..."<<endl;
	cout<<endl;
	CPA.TrackingClusters(R);
	cout<<"Done"<<endl;
	cout<<endl;

	cout<< "Visualizing global topic weights..."<<endl;

	if(N<10001 && Iterations < 100)  // ho salvato It beta e It partizioni, disegno la partizione in ogni caso
		CPA.VisualizeBeta(R,BestClustering);

	else if(N<10001 && Iterations >= 100 && Iterations < 10000){ // ho salvato 100 Beta e It partizioni, devo capire se disegnare la partizione
		if(Iterations - BestClustering <= 100)
			CPA.VisualizeBeta(R,BestClustering - Iterations + 100);   // ho trovato il miglior clustering nelle ultime 100 iterazioni
		else
			CPA.VisualizeBeta(R,std::numeric_limits<unsigned long>::max());
	}

	else if(N<10001 && Iterations >= 10000){ // ho salvato 100 Beta e 10000 partizioni, devo capire se disegnare la partizione
		if(10000 - BestClustering <= 100)
			CPA.VisualizeBeta(R,BestClustering - 10000 + 100);   // ho trovato il miglior clustering nelle ultime 100 iterazioni
		else
			CPA.VisualizeBeta(R,std::numeric_limits<unsigned long>::max());
	}

	else
        CPA.VisualizeBeta(R,std::numeric_limits<unsigned long>::max());

	cout<<"Done"<<endl;
	cout<<endl;

    cout<<"Posterior analysis of Numbers of Cluster"<<endl;
	cout<<endl;

	while(Try == 'y'){

		cout<<"Burnin (choose a value <= "<<Kdim<<"): "<<endl;
		cin>>Burnin;
		while(Burnin > Kdim){
			cout<<"Burnin (choose a value <= "<<Kdim<<"): "<<endl;
			cin>>Burnin;
		}
		cout<<"Thinning: (must be > 0)"<<endl;
		cin>>Thinning;
		while(Thinning <=0){
		cout<<"Thinning: (must be > 0)"<<endl;
		cin>>Thinning;
		}

		cout<<"Doing the posterior analysis"<<endl;
		CPA.KPosteriorAnalysis(R,Burnin, Thinning);
		cout<<"Look in R_results directory if your analysis is good or you would like to change something"<<endl;
		cout<<"..."<<endl;
		cout<<"Would you like to change your posterior analysis?  press y[yes] or n[no]  "<<endl;
		cin>>Try;

		while(Try !='y' && Try != 'n'){
			cout<<"press y [yes] or n [no] "<<endl;
			cin>>Try;
		}

	}

	cout<<endl;

	if(Prior == 'y'){
		Try = 'y';

		std::cout<<"Posterior Analysis of Gamma and Alpha"<<std::endl;

		while(Try == 'y'){

			if(AlphaTry == 'y' ){
				cout<<"ANALYSIS ON ALPHA"<<endl;
				cout<<"Burnin (choose a value <= "<<Alphadim<<"): "<<endl;
				cin>>AlphaBurnin;

				while(AlphaBurnin > Alphadim){
					cout<<"Burnin (choose a value <= "<<Alphadim<<"): "<<endl;
					cin>>AlphaBurnin;
				}
				cout<<"Thinning: (must be > 0)"<<endl;
				cin>>AlphaThinning;

				while(AlphaThinning <=0){
					cout<<"Thinning: (must be > 0)"<<endl;
					cin>>AlphaThinning;
				}
			}
			if(GammaTry == 'y' ){
				cout<<"ANALYSIS ON GAMMA"<<endl;
				cout<<"Burnin (choose a value <= "<<Gammadim<<"): "<<endl;
				cin>>GammaBurnin;

				while(GammaBurnin > Gammadim){
					cout<<"Burnin (choose a value <= "<<Gammadim<<"): "<<endl;
					cin>>GammaBurnin;
				}
				cout<<"Thinning: (must be > 0)"<<endl;
				cin>>GammaThinning;

				while(GammaThinning <=0){
					cout<<"Thinning: (must be > 0)"<<endl;
					cin>>GammaThinning;
				}

			}

			cout<<"Doing the Alpha Gamma posterior analysis..."<<endl;
			CPA.AGPosteriorAnalysis(R,AlphaBurnin, AlphaThinning, GammaBurnin, GammaThinning, AlphaTry, GammaTry);


			cout<<"Look in R_results directory if your analysis is good or you would like to change something"<<endl;
			cout<<"..."<<endl;
			cout<<"Would you like to change your posterior analysis?  press y[yes] or n[no]  "<<endl;
			cin>>Try;

			while(Try !='y' && Try != 'n'){
				cout<<"press y [yes] or n [no] "<<endl;
				cin>>Try;
			}

			if(Try == 'y'){
				cout<<"Would you like to change Alpha? press y [yes] or n [no]  "<<endl;
				cin>>AlphaTry;

				while(AlphaTry !='y' && AlphaTry != 'n'){
					cout<<"press y [yes] or n [no] "<<endl;
					cin>>AlphaTry;
				}

				cout<<"Would you like to change Gamma? press y [yes] or n [no]  "<<endl;
				cin>>GammaTry;

				while(GammaTry !='y' && GammaTry != 'n'){
					cout<<"press y [yes] or n [no] "<<endl;
					cin>>GammaTry;
				}
			}


		}

	}

	cout<<endl;

	if(is_LPML){
        cout<<"Computing LPML..."<<endl;
			CPA.LPML(Iterations-burnin);
	}
	cout<<"Done"<<endl;
	cout<<endl;
	*/

	return 0;
}
