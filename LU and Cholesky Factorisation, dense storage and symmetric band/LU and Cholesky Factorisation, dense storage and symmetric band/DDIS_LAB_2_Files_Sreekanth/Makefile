# compilation flag for optimized (Release mode)
CXXFLAGS = -O3 -std=c++11
#compiler
CXX = g++
# compilation flag for debug mode. To use comment the above line and uncomment the line below
#CXXFLAGS = -g -std=c++11

# BLASLIB is the flag to link with your BLAS implementation. below is the default: the system bls lib. 
#BLASLIB = -lblas -llapack
# To use an other blas complete and uncomment the line below, set for openblas
# You can get the openblas source from here : 
# www.openblas.net, get the tar.gz
#  copy it in your directory 
#  extract it : tar -zxvf OpenBLAS-0.2.19.tar.gz
#  compile it : cd OpenBLAS-0.2.19
#               make
#  When it is done you should have the file libopenblas.a, the openblas library
#  export OPENBLAS_NUM_THREADS=2
# replace OpenBLAS-0.2.19/libopenblas.a to your actuall path of your openblaslibrary.


#BLASLIB =  libopenblas.a -lpthread -lgfortran
BLASLIB = /home/srbakkache2022/Desktop/M2/LAB_1/openblas/OpenBLAS-0.3.24/libopenblas.a  -lpthread -lgfortran

all: lab2_part1 lab2_part2

# this are the rule to compile the library file that Igave you (mmio.cc and dmatrix_denseCM.cc)
dmatrix_denseCM.o : Makefile dmatrix_denseCM.cc mmio.h dmatrix_denseCM.h
	$(CXX) $(CXXFLAGS) -c dmatrix_denseCM.cc

mmio.o :  Makefile mmio.cc mmio.h
	$(CXX) $(CXXFLAGS) -c mmio.cc	

dsquarematrix_symband.o : Makefile dsquarematrix_symband.cc mmio.h dsquarematrix_symband.h
	$(CXX) $(CXXFLAGS) -c dsquarematrix_symband.cc

dLU_denseCM.o : Makefile dLU_denseCM.cc mmio.h dmatrix_denseCM.h dLU_denseCM.h 
	$(CXX) $(CXXFLAGS) -c dLU_denseCM.cc

dCholesky_band.o: Makefile dCholesky_band.cc mmio.h dCholesky_band.h dmatrix_denseCM.h dsquarematrix_symband.h blas_lapack_def.h 
	$(CXX) $(CXXFLAGS) -c dCholesky_band.cc


lab2_part1.o : Makefile lab2_part1.cc dmatrix_denseCM.h dLU_denseCM.h
	$(CXX) $(CXXFLAGS) -c lab2_part1.cc

lab2_part2.o : Makefile lab2_part2.cc dmatrix_denseCM.h dLU_denseCM.h dsquarematrix_symband.h dCholesky_band.h
	$(CXX) $(CXXFLAGS) -c lab2_part2.cc


lab2_part1 : Makefile lab2_part1.o dmatrix_denseCM.o mmio.o dLU_denseCM.o
	$(CXX) $(CXXFLAGS) lab2_part1.o dmatrix_denseCM.o mmio.o dLU_denseCM.o $(BLASLIB) -o lab2_part1

lab2_part2 : Makefile lab2_part2.o dmatrix_denseCM.o dsquarematrix_symband.o dCholesky_band.o mmio.o dLU_denseCM.o
	$(CXX) $(CXXFLAGS) lab2_part2.o dmatrix_denseCM.o mmio.o dLU_denseCM.o dsquarematrix_symband.o dCholesky_band.o $(BLASLIB) -o lab2_part2

clean : 
	rm -rf *.o *~
	rm lab2_part1 lab2_part2



