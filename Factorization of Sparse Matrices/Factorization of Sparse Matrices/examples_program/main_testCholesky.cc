#include "dLU_denseCM.h"
#include "dLLT_denseCM.h"
#include <chrono>
#include "mesh.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <random>
#include <functional>

dmatrix_denseCM generateRandom(int m, int n){
  dmatrix_denseCM A(m,n);
  std::random_device rd;
std::mt19937 gen(rd());  // this is the generator
 std::uniform_real_distribution<> dis(0., 1.); // this is to have an uniform distribution between 0. and 1.
 std::generate(A.begin(), A.end(), std::bind(dis,gen));
 return A;
}

int main(){
  
  
  std::string inputfile ="data/square4.msh";
  std::ifstream meshin(inputfile);
  if (!meshin.is_open()) {std::cout << "Can't open mesh file " << inputfile << std::endl;   throw;}
  std::cout << "Reading Mesh from file " << inputfile << std::endl; 
  mesh me;
  meshin >> me;
  std::cout << me.nodes.size() << " "<< me.triangles.size() << std::endl;
  fem fem1(me);
  dmatrix_denseCM A = fem1.generateMassMatrix_denseCM();
  
  /*
  std::string inputfilename ="data_band.mat";
  std::ifstream inputfile(inputfilename);
  dmatrix_denseCM A;
  inputfile >> A;
  //std::cout << A << std::endl;
  */
  int m = A.getNbLines();
  dmatrix_denseCM b = generateRandom(m,1);

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  dLLT_denseCM LLT_A(A);
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "time LLT " <<  (elapsed_seconds).count() << std::endl;
  
  //  std::cout << LLT_A.getL() << std::endl;
  dmatrix_denseCM LLT = LLT_A.getL() * transpose(LLT_A.getL());
  //std::cout << LLT << std::endl;
  
  dmatrix_denseCM x = LLT_A.solve(b);
  std::cout << norm2(A*x - b) << std::endl;
  
 
  
}
