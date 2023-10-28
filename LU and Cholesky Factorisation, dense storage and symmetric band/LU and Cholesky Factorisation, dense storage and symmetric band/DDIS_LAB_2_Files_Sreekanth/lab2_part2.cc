#include <algorithm>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include<functional>
#include<math.h>

#include "dmatrix_denseCM.h"
#include "dLU_denseCM.h"
#include "dCholesky_band.h"

int main(){
  std::cout << "m, elapsed_time " << std::endl;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.1, 1.);
  std::chrono::time_point<std::chrono::system_clock> start, end;
  
  /// If verbose == true, all the matrices will be outputed to the screen for verification. Set it to false on large matrix !
  const bool verbose = false;


  dmatrix_denseCM A; // A object of class dmatrix_denseCM is created.
  // data_band.mat contain a band matrix that you can use for your test
  //std::string filename("data_band.mat");
  //bcsstk14.mtx contain a relatively large positive definite band band matrix  
  //std::string filename("bcsstk14.mtx");
  //std::string filename("bcsstk15.mtx"); // a larger matrix.
  //std::string filename("bcsstk16.mtx");
  //std::string filename("bcsstk17.mtx");
  //std::string filename("bcsstk18.mtx"); // a larger matrix.
  std::string filename("bcsstk27.mtx");

  
  std::ifstream in(filename.c_str());
 
  if(!in.is_open() ) std::cout << "Can't Open File " << filename  << std::endl;
  in >> A; //the data from the file is read into A
  start = std::chrono::system_clock::now();
  dLU_denseCM LU(A, factorpolicy(factorpolicy::Lapack)); // using lapack library the LU decomposition is computed, to compare our computed solution.
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "Time for factorising mat from  " << filename  << ": " << elapsed_seconds.count() << "s" << std::endl;
  
  dmatrix_denseCM B(A.getNbLines(), 1); 
  std::generate(B.begin(), B.end(), std::bind(dis,gen) );
  /// This should give you the correct answer, since it used the factorisation algorithm that I gave you
  dmatrix_denseCM X = LU.solve(B);
  std::cout << " Error in solving the linear system " << norm2(A*X-B) << std::endl;
 
  // This get the up and down bandwith of A. You need to code this functions ...
  int lbu = computeBandwidthUp(A);
  int lbd = computeBandwidthDown(A);
  assert (lbu == lbd); // A is supposed symmmetric, so that lbu should be equal to lbd
  std::cout << "lbu " << lbu << std::endl;
  std::cout << "lbd " << lbd << std::endl;
  // This copy in the symmetric band storage Ab.
  dsquarematrix_symband Ab(A, lbu);
  if(verbose) {
   std::cout << "Checking Band Matrix Storage :" << std::endl;
   std::cout << "A in dense storage" << std::endl;
   std::cout << A << std::endl;
   std::cout << "Ab : A in compressed symetric band storage (storing the upper part of the matrix.)" << std::endl;
   std::cout << Ab << std::endl;
  }

  start = std::chrono::system_clock::now();
  // This compute the Cholesky Factorisation of the symmetric positive definite band matrix stored in Ab, and store it in cholAB.
  dCholesky_band cholAb(Ab);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << " Cholesky band factorisation Done,  Time " <<  elapsed_seconds.count()   << " seconds " << std::endl;
  // This solve the linear system AX2 = B using the cholesky band_ factorisation.
  dmatrix_denseCM X2 =  cholAb.solve(B);
  /// If your implementation of what I ask for is correct, this should give you something very close to 0 : We are comparing the soltion using the full and the band storage.
  std::cout << " Error in the factorisation " << norm2(X-X2) << std::endl;
}

