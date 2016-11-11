
IFR	:=	useR


SRC :=				rngstream.cpp omprng.cpp
OBJ :=				../src/include/rngstream.o ../src/include/omprng.o
OBJFILE :=			rngstream.o omprng.o
MAIN :=				main.cpp
MAINDTS :=			mainDataset.cpp
EXECDTS :=			dataset
EXEC :=				main
DTS :=				Dataset.txt Variables.txt Vocabulary.txt

#path to find include file to compile them and cancel the obj file
SRCPATH :=		cd ../src/include;
#path to find output of the algorithm. They are file.txt
CPPOUTPATH :=	cd cpp_results;

#flags that is common
OTHERFLAGS :=		-Wall -std=c++11 -fopenmp


ifeq ($(IFR), useR)

R_HOME := 		$(shell R RHOME)

## include headers and libraries for R
RCPPFLAGS :=	$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS :=	$(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := 	$(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK :=	$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

## if you need to set an rpath to R itself, also uncomment
#RRPATH :=		-Wl,-rpath,$(R_HOME)/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := 	$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 	$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


## include headers and libraries for RInside embedding classes
RINSIDEINCL :=	$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS :=	$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## compiler etc settings used in default make rules
CXX :=	$(shell $(R_HOME)/bin/R CMD config CXX)
CPPFLAGS :=	$(shell $(R_HOME)/bin/R CMD config CPPFLAGS)
CXXFLAGS :=	$(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS)
LDLIBS :=	$(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)
LDFLAGS :=	-L${mkBlasLib} -lopenblas
#OTHERFLAGS :=		-Wall -std=c++11

## path for include files
#il makefile si trova nella cartella del main, qui salto ad un livello superiore
OUTPUTPATH :=	cd R_results;


#compile only main program
build_main:
			rm -f *~ $(EXEC)
			$(OUTPUTPATH) rm -f *
			$(CPPOUTPATH) rm -f *
			$(SRCPATH) rm -f $(OBJFILE)
			$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
			$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OTHERFLAGS) $(MAIN) $(LDLIBS) $(LDFLAGS) -o $(EXEC) $(OBJ) -DPA

#compile only dataset
build_dataset:
					rm -f *~ $(EXECDTS) $(DTS)
					$(OUTPUTPATH) rm -f *
					$(CPPOUTPATH) rm -f *
					$(SRCPATH) rm -f $(OBJFILE)
					$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
					$(CXX) $(OTHERFLAGS) $(MAINDTS) -o $(EXECDTS) $(OBJ)


#compile and run dataset and main
run_all:
			rm -f *~ $(EXEC) $(EXECDTS) $(DTS)
			$(OUTPUTPATH) rm -f *
			$(CPPOUTPATH) rm -f *
			$(SRCPATH) rm -f $(OBJFILE)
			$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
			$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OTHERFLAGS) $(MAIN) $(LDLIBS) $(LDFLAGS) -o $(EXEC) $(OBJ) -DPA
			$(CXX) $(OTHERFLAGS) $(MAINDTS) -o $(EXECDTS) $(OBJ)
			./$(EXECDTS) -f; ./$(EXEC) -f

#compile and run only main program
run_main:
			rm -f *~ $(EXEC)
			$(OUTPUTPATH) rm -f *
			$(CPPOUTPATH) rm -f *
			$(SRCPATH) rm -f $(OBJFILE)
			$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
			$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OTHERFLAGS) $(MAIN) $(LDLIBS) $(LDFLAGS) -o $(EXEC) $(OBJ) -DPA
			./$(EXEC) -f

#compile and run only dataset
run_dataset:
					rm -f *~ $(EXECDTS) $(DTS)
					$(OUTPUTPATH) rm -f *
					$(CPPOUTPATH) rm -f *
					$(SRCPATH) rm -f $(OBJFILE)
					$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
					$(CXX) $(OTHERFLAGS) $(MAINDTS) -o $(EXECDTS) $(OBJ)
					./$(EXECDTS) -f



#delete executable
clean:
			rm -f *~ $(EXEC) $(EXECDTS)

#delete executable, objects, dataset and results
distclean:
			rm -f *~ $(EXEC) $(EXECDTS) $(DTS)
			$(OUTPUTPATH) rm -f *
			$(CPPOUTPATH) rm -f *
			$(SRCPATH) rm -f $(OBJFILE)

else

# make IFR="notuseR"

#compile only main program
build_main:
			rm -f *~ $(EXEC)
			$(CPPOUTPATH) rm -f *
			$(SRCPATH) rm -f $(OBJFILE)
			$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
			$(CXX) $(OTHERFLAGS) $(MAIN) -o $(EXEC) $(OBJ)

#compile and run dataset and main
run_all:
			rm -f *~ $(EXEC) $(EXECDTS) $(DTS)
			$(CPPOUTPATH) rm -f *
			$(SRCPATH) rm -f $(OBJFILE)
			$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
			$(CXX) $(OTHERFLAGS) $(MAIN) -o $(EXEC) $(OBJ)
			$(CXX) $(OTHERFLAGS) $(MAINDTS) -o $(EXECDTS) $(OBJ)
			./$(EXECDTS) -f;./$(EXEC) -f
#run only main program
run_main:
			rm -f *~ $(EXEC)
			$(CPPOUTPATH) rm -f *
			$(SRCPATH) rm -f $(OBJFILE)
			$(SRCPATH) $(CXX) $(OTHERFLAGS) -c $(SRC)
			$(CXX) $(OTHERFLAGS) $(MAIN) -o $(EXEC) $(OBJ)
			./$(EXEC) -f


endif
