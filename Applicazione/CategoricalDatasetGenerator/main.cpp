#include "../src/include/HDP_MCMC.hpp"
#include "../src/include/Document.hpp"
#include "../src/include/Model.hpp"
#include "../src/include/Cluster.hpp"
#include "../src/include/Type.hpp"
#include "../src/include/Functions.hpp"

#include <iostream>
#include <string>

#ifdef PA
    #include "../src/include/PosteriorAnalysis.hpp"
    #include <RInside.h>
#endif

#include "../src/include/GetPot"


int main(int argc, char** argv){

    using std::cout;
    using std::endl;
    using std::cin;
    using std::cerr;

    GetPot command_line (argc,argv);
    const std::string DataFile = command_line.follow("DataFile.txt",1,"-f");

    GetPot file (DataFile.c_str());
    //Alpha
    bool is_alpha = static_cast<bool>(file("is_alpha", 1));
    bool is_alpha_prior = static_cast<bool>(file("is_alpha_prior",0));
    double Alpha = file("a", 1.0);
    double AA = file("aa", 5.0);
    double AB = file("ab", 0.1);
    //Gamma
    bool is_gamma = static_cast<bool>(file("is_gamma",0));
    bool is_gamma_prior = static_cast<bool>(file("is_gamma_prior",0));
    double Gamma = file("g", 1.0);
    double GA = file("ga", 5.0);
    double GB = file("gb",0.1);
    //Lambda
    bool is_lambda = static_cast<bool>(file("is_lambda",1));
    double l = file("l", 1.5);
    unsigned int Dim = file("Dim",3);
    vector<double> Lambda(Dim,l);
	//Seed
	bool is_seed = static_cast<bool> (file ("is_seed", 0));
    const unsigned long Seed = file("seed", 0);

	//K iniziale
	bool is_K_init= static_cast<bool> (file ("is_K_init", 1));
	unsigned int K_init = file("K_init",10);

	//Iteration
	const unsigned long Iterations = file("Iter",10000);
	const unsigned long SubIterations = file("SubIter",10);

	bool is_LPML = static_cast<bool> (file ("is_LPML", 0));
	const unsigned long burnin = file("burnin",0);

    const unsigned int DIM=1;
    HDP_MCMC <CategoricalModel,CategoricalDocument,DIM> sm;

	cout<<endl;
    cout<<"HDP_MCMC"<<endl;
	cout<<endl;

    if(is_alpha){
	    cout<<"Fixed Alpha"<<endl;
		cout<<"Alpha = "<<Alpha<<endl;
        sm.SetAlphaFixed(Alpha);
    }
    else{
        if(is_alpha_prior){
			cout<<"Prior on Alpha"<<endl;
			cout<<"Hyperparameters: a = "<<AA<<" e b = "<<AB<<endl;
            sm.SetAlphaPrior(AA,AB);
        }
    }

	cout<<endl;

    if(is_gamma){
        cout<<"Fixed Gamma"<<endl;
		cout<<"Gamma = "<<Gamma<<endl;
        sm.SetGammaFixed(Gamma);
    }
    else{
        if(is_gamma_prior){
			cout<<"Prior on Gamma"<<endl;
			cout<<"Hyperparameters: a = "<<GA<<" e b = "<<GB<<endl;
            sm.SetGammaPrior(GA,GB);
    }
    }

	cout<<endl;

    if(is_lambda){
        cout<<"Dim Lambda = "<<Dim<<endl;
        sm.SetLambdaInfo(Lambda);
    }

	cout<<endl;

	if(is_seed){
        cout<<"Seed = "<<Seed<<endl;
        sm.SetSeed(Seed);
    }

	cout<<endl;

	if(is_K_init){
	  if(K_init>Dim){
		cout<<"You can't have a number of clusters higher than the number of distinct words"<<endl;
		cout<<"Choose a number of cluster K<"<<Dim<<": "<<endl;
		cin>>K_init;
		while(K_init>Dim){
			cout<<"Choose a number of clusters K<"<<Dim<<": "<<endl;
			cin>>K_init;
		}
	  }
	  else
	  	cout<< "Initial K = "<<K_init<<endl;
	  sm.SetK_init(K_init);
	}
	else{
		cout<<"Initial K is generated randomly between 1 and 10"<<endl;
	}

	cout<<endl;

	if(is_LPML){
		cout<<"You are monitoring the model with LPML"<<endl;
		sm.Check_Model(burnin);
	}

	cout<<endl;

    cout<<"Loading Dataset..."<<endl;
    sm.SetDataset("Dataset.txt", "Variables.txt");
	cout<<"Done"<<endl;
	cout<<endl;

    cout<<"Starting..."<<endl;

    double start = omp_get_wtime();
    sm.Algorithm(Iterations,SubIterations);
    double stop = omp_get_wtime();
	cout<<"Done"<<endl;
	cout<<endl;

    cout<< "Execution time : "<< stop-start <<" sec "<<endl;
	cout<<endl;


#ifdef PA

    cout<<"Starting Posterior Analysis..."<<endl;
    cout<<endl;

    // variabili utilizzate nell'analisi a posteriori
    const std::string wd="setwd(\"~/APSC_Parisi_Perego/CategoricalDatasetGenerator/R_results\")";		// imposta la wd per R, la passo ad ogni funzione che usa R
	unsigned int W;   //  nr parole distinte
    unsigned int D;   // nr documenti
	unsigned int N;   // numerosità dataset
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

    // Occorre un oggetto RInside che viene passato per referenza alla funzione
    cout <<"Opening R environment... "<<endl;
    RInside R(argc, argv);
    cout <<"Done. " <<endl;
	cout<<endl;

	W = sm.ViewW();
    D = sm.ViewD();
	N = sm.ViewN();

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
	
	unsigned long BestClustering = 0;
	
	if(N<10001){
		cout<<"Loading labels..."<<endl;
		if(Iterations >= 10000)
			CPA.LoadLabels(10000);
		else
			CPA.LoadLabels(Iterations);
		cout<<"Done"<<endl;
		cout<<endl;

		cout<<"Finding least square clustering..."<<endl;
		cout<<endl;
		BestClustering = CPA.LeastSquareClustering();
		cout<<"Done"<<endl;
		cout<<endl;
	}

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

	cout<<"Loading all K..."<<endl;
	Kdim = CPA.SetAllK();
	cout<<"Done."<<endl;
	cout<<endl;

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

	cout<<"Loading vocabulary..."<<endl;
	CPA.SetVocabulary();
	cout<<"Done"<<endl;
	cout<<endl;

	cout<< "Loading theta..."<<endl;
    CPA.SetTheta(Iterations);
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
				cout<<"GAMMA ON ALPHA"<<endl;
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


#endif

	return 0;
}
