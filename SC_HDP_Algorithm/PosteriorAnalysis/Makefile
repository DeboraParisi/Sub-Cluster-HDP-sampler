R_HOME := 		$(shell R RHOME)


SRC :=				rngstream.cpp omprng.cpp
OBJ :=				../src/include/rngstream.o ../src/include/omprng.o
OBJFILE :=			rngstream.o omprng.o
#MAIN :=				PosteriorAnalysis.cpp
#MAINDTS :=			Generator.cpp
#EXECDTS :=			dataset
EXEC :=				analysis
#this variable need to set the directory where are the files to analays
DIRECTORY ?= PosteriorAnalysys

## include headers and libraries for R
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := 		$(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK := 		$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

## if you need to set an rpath to R itself, also uncomment
#RRPATH :=		-Wl,-rpath,$(R_HOME)/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


## include headers and libraries for RInside embedding classes
RINSIDEINCL := 		$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := 		$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## compiler etc settings used in default make rules
CXX := 			$(shell $(R_HOME)/bin/R CMD config CXX)
CPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config CPPFLAGS)
CXXFLAGS := 		$(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS)
LDLIBS := 		$(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)
LDFLAGS := 		-L${mkBlasLib} -lopenblas
OTHERFLAGS :=		-Wall -std=c++11 -fopenmp
#OTHERFLAGS :=		-Wall -std=c++11

## path for include files
#il makefile si trova nella cartella del main, qui salto ad un livello superiore
SRCPATH :=		cd ../src/include;
OUTPUTPATH :=		cd R_results;
CPPOUTPATH :=		cd cpp_results;

ifeq ($(DIRECTORY), PosteriorAnalysys)

all:
	@echo "Error: You haven't set a Directory in which you can find the files to analyse"
	exit
main:
	@echo "Error: You haven't set a Directory in which you can find the files to analyse"
	exit

else

#compile
build_categorical:
			rm -f *~ $(EXEC)
			$(OUTPUTPATH) rm -f *
			$(CPPOUTPATH) rm -f *
			cp ../$(DIRECTORY)/cpp_results/*.txt ./cpp_results
		    cp ../$(DIRECTORY)/cpp_results/*.bin ./cpp_results
			cp ../$(DIRECTORY)/Vocabulary.txt .
			cp ../$(DIRECTORY)/Variables.txt .
			cp ../$(DIRECTORY)/DataFile.txt .
			#cp ../$(DIRECTORY)/NewTitolo.txt .
			$(SRCPATH) rm -f $(OBJFILE)
			$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
			$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OTHERFLAGS) Categorical.cpp $(LDLIBS) $(LDFLAGS) -o $(EXEC) $(OBJ)



#compile and run
run_categorical:
			rm -f *~ $(EXEC)
			$(OUTPUTPATH) rm -f *
			$(CPPOUTPATH) rm -f *
			cp ../$(DIRECTORY)/cpp_results/* ./cpp_results
			cp ../$(DIRECTORY)/cpp_results/*.bin ./cpp_results
			cp ../$(DIRECTORY)/Vocabulary.txt .
			cp ../$(DIRECTORY)/Variables.txt .
			cp ../$(DIRECTORY)/DataFile.txt .
			$(SRCPATH) rm -f $(OBJFILE)
			$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
			$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OTHERFLAGS) Categorical.cpp $(LDLIBS) $(LDFLAGS) -o $(EXEC) $(OBJ)
			./$(EXEC)

endif


#delete executable and analysis
distclean:
			rm -f *~ $(EXEC)
			rm -f *.txt
			rm -f *.bin
			$(OUTPUTPATH) rm -f *
			$(CPPOUTPATH) rm -f *

