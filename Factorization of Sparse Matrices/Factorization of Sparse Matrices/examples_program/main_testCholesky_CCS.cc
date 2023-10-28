#include "dLU_denseCM.h"
#include "dLLT_denseCM.h"
#include "dmatrix_CCS.h"
#include "dLLT_CCS.h"

#include <chrono>
#include "mesh.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <random>

#include "mat2gnuplot.h"
#include "dmatrix_denseCM.h"



dmatrix_CCS generateFromCommandLineArguments(int argc, char *argv[]){
  std::string inputfile = "data/small_sym_sparse.mtx";
  bool printmat  = false;
  
  if (argc < 2){
    std::cout << "Usage : " << std::endl;
    std::cout << argv[0] << " input_filename" <<  " 0 " << std::endl;
    std::cout << "or :" << std::endl;
    std::cout << argv[0] << " input_filename" <<  " 1 " << std::endl;
    std::cout << std::endl<<  " Where input_filename is the input file, which can be a matrix format file (.mtx) or a gmsh mesh file (.msh ) " << std::endl;
    std::cout << " - If a mesh file is given, the analysis will be done on the classical finite element mass matrix build from the mesh " << std::endl;
    std::cout << " - If a mtx file is given, the analysis will be done directly on the matrix discribed by this file" << std::endl;
    std:: cout << " The second argument (0 or 1), is used to tell the code if you want to print the matrix shape on the screen (1), or not" << std::endl;
    std::cout << "Trying to run with default value : " << inputfile << " " << printmat << std::endl;
  }
  else{
    inputfile = argv[1];
    if (argc > 2)
      printmat = atoi(argv[2]);
  }
  bool from_mesh = (inputfile.npos != inputfile.find(".msh" ) );
  bool from_mtx = (inputfile.npos != inputfile.find(".mtx" ) );
  
  dmatrix_CCS A;
  if (from_mtx){
    std::ifstream matfile(inputfile);
    if (!matfile.is_open()){ std::cout << "Can't open matrix file " << inputfile << std::endl; throw;}
    std::cout << "Reading Matrix from file " << inputfile << std::endl; 
    matfile >> A;
  }
  if (from_mesh){
    std::ifstream meshin(inputfile);
    if (!meshin.is_open()) {std::cout << "Can't open mesh file " << inputfile << std::endl;   throw;}
    std::cout << "Reading Mesh from file " << inputfile << std::endl; 
    mesh m;
    meshin >> m;
    std::cout << m.nodes.size() << " "<< m.triangles.size() << std::endl;
    fem fem1(m);
    A = fem1.generateMassMatrix_CCS();
  }
  
  if (printmat) std::cout << "A : \n" <<  convertToDenseCM(A) << std::endl;
  std::cout << "matrix size : m = " << A.m<< " n = " << A.n << " nnz =" << A.nnz << std::endl;
  return A;
}

int main(int argc, char *argv[]){
  // generate a matrix from command line argument. Since we are using Cholesky decomposition here, A should be symetric postiv definite
  dmatrix_CCS A = generateFromCommandLineArguments(argc, argv);
  
  // Variable to store the start and end of a sequence of code for timing purpose
  std::chrono::time_point<std::chrono::system_clock> start, end;
  
  int m = A.getNbLines();
  // B is the Right hand side matrix to solve for. here a random matrix, with m line and 3 column (we want to solve the system fo 3 right and side )
  int nrhs = 3; 
  dmatrix_denseCM B = generateRandom(m,nrhs, 0., 1.);
  // Xref will store the reference result, the solution of AX =B;
  dmatrix_denseCM Xref = (m, nrhs, 0.);
  std::cout << "Trying to solve the original system " << std::endl;
  // Compute the Cholesky Factorisation of A
  start = std::chrono::system_clock::now();
  dLLT_CCS LLT_A(A);
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << " Time Cholesky sparse " <<  (elapsed_seconds).count() << std::endl;
  // Solve AX= B using the factorisation of A.
  Xref = LLT_A.solve(B);
  std::cout << " NNZ in L " << LLT_A.getL().nnz << " " << "LB of  L : " <<  lbdown(LLT_A.getL()) << std::endl;
  // Checking result. If solve is correct, result should be of Norm close to Zero
  std::cout << " Checking results : norm2(A*X -B) " << norm2(convertToDenseCM(A)*Xref - B) << std::endl;
  
  // Now we try again, after permuting the system 
  std::vector< int > neworder = cuthillmckee( buildSymGraph(A) );
  { 
    std::cout << "Trying to solve the permuted system using cuthillmckee ordering " << std::endl;
    /// permuting the system  :  pA = PAPT, pB = PB
    dmatrix_CCS     pAp = sympermut(A, neworder);
    dmatrix_denseCM pB =  permutlines(B, neworder);
    start = std::chrono::system_clock::now();
    
    dLLT_CCS LLT_pAp(pAp);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << " NNZ in L " << LLT_pAp.getL().nnz << " " << "LB of  L : " <<  lbdown(LLT_pAp.getL()) << std::endl;
    std::cout << " time Cholesky Sparse  " <<  (elapsed_seconds).count() << std::endl;
    dmatrix_denseCM pX = LLT_pAp.solve(pB);
    
    
    
    dmatrix_denseCM AT = convertToDenseCM(LLT_pAp.getLT());
    dmatrix_denseCM A =  convertToDenseCM(LLT_pAp.getL());
    dmatrix_denseCM B = A + AT;
    //spy_ge(B, "CUTHUILL_CUBE2.plt" );
    //spy_ge(B, "CUTHUILL_CUBE.plt" );
    //spy_ge(B, "CUTHUILL_sq.plt" );
    //spy_ge(B, "CUTHUILL_sq2.plt" );
    //spy_ge(B, "CUTHUILL_sq3.plt" );
    //spy_ge(B, "CUTHUILL_square4.plt" );
    //spy_ge(B, "CUTHUILL_b14.plt" );
    spy_ge(B, "CUTHUILL_b15.plt" );
    
    
    
    
    
    
    
    
    
    
    
    
    
    dmatrix_denseCM X  = ipermutlines(pX, neworder);
    std::cout << " check : error Norm2(X-Xref) " << norm2(X-Xref) << std::endl;
  }
  { 
    std::cout << "Trying to solve the permuted system using reverse cuthillmckee ordering " << std::endl;
    std::vector< int > rorder(neworder.size());
    std::copy(neworder.rbegin(), neworder.rend(), rorder.begin());
    /// permuting the system  :  pA = PAPT, pB = PB
    dmatrix_CCS     pAp = sympermut(A, rorder);
    dmatrix_denseCM pB =  permutlines(B, rorder);
    start = std::chrono::system_clock::now();
    dLLT_CCS LLT_pAp(pAp);
    end = std::chrono::system_clock::now();
    
    
    

    dmatrix_denseCM AT = convertToDenseCM(LLT_pAp.getLT());
    dmatrix_denseCM A =  convertToDenseCM(LLT_pAp.getL());
    dmatrix_denseCM B = A + AT;
    //spy_ge(B, "ReverseCUTHUILL_CUBE2.plt" );
    //spy_ge(B, "ReverseCUTHUILL_CUBE.plt" );
    //spy_ge(B, "ReverseCUTHUILL_square.plt" );
    //spy_ge(B, "ReverseCUTHUILL_square2.plt" );
    //spy_ge(B, "ReverseCUTHUILL_square3.plt" );
    //spy_ge(B, "ReverseCUTHUILL_square4.plt" );
    
    //spy_ge(B, "ReverseCUTHUILL_b14.plt" );
    spy_ge(B, "ReverseCUTHUILL_b15.plt" );
    
    
    
    
    
    
    
    elapsed_seconds = end-start;
    std::cout << " NNZ in L " << LLT_pAp.getL().nnz << " " << "LB of  L : " <<  lbdown(LLT_pAp.getL()) << std::endl;
    std::cout << " time Cholesky Sparse  " <<  (elapsed_seconds).count() << std::endl;
    dmatrix_denseCM pX = LLT_pAp.solve(pB);
    dmatrix_denseCM X  = ipermutlines(pX, rorder);
    std::cout << " check : error Norm2(X-Xref) " << norm2(X-Xref) << std::endl;
  }
  
}
