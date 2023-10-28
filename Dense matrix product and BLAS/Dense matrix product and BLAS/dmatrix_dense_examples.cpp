// This program contain examples to show you some usefull syntaxes and construct. This is here to help you build your hown.

// A program star by including the definitiions of what you need.

// You need to include the following file to use the dmatrix_dense_CM class in your programm.
#include "dmatrix_denseCM.h" 
// this is included so that you can use std::cout, std::endl; ... to print stuff on the screen
#include <iostream>
// This is to be able to open files.
#include <fstream>
// This is to generate random numbers
#include <random>
// This contain stl algorithm like generate, copy etc ..
#include <algorithm>
#include<functional>

int main(){
  // Create a Matrix A of size (4*5), filled with 0.; 
  dmatrix_denseCM A(4,5);
  // print the matrix on screen :
  std::cout << A << std::endl;
  // print the first  entry of the matrix  (0,0):
  std::cout << A(0,0) << std::endl;
  // print the last entry (3,4):
  std::cout << A(3,4) << std::endl;
  // set entry (3,4) to 25;
  A(3,4) = 25;
  // print the matrix on screen :
  std::cout << A << std::endl;
  // get the number of lines and columns of the matrix :
  std::cout << A.getNbLines() << " " << A.getNbColumns() << std::endl;
  // get the norm2 of the matrix (sqrt(sum_ij (A(ij)^2)))
  std::cout << A.norm2() << std::endl;
  // Create 2 empty matrix B, C:
  dmatrix_denseCM B, C;
  // Read  matrix B from the screen
  std::cin >> B;
  // print the B matrix on screen :
  std::cout << "You just entered the following Matrix (B)" << std::endl;
  std::cout << B<< std::endl;
  // Read  matrix C from the screen
  std::cin >> C;
  // print the B matrix on screen :
  std::cout << "You just entered the following Matrix (C)" << std::endl;
  std::cout << B<< std::endl;
  // Try to add B and C into the new matrix D
  dmatrix_denseCM D = B+C;
  //print the D matrix;
  std::cout << "D = B+C" << std::endl << D << std::endl;
  // substract C from D
  D -= C;
  //print D
  std::cout << D << std::endl;
  // multiply each terms of D by 2;
  D *=2.;
  //print D
  std::cout << D << std::endl;
  // divide each terms of D by 2;
  D /=2.;
  //print D
  std::cout << D << std::endl;
  // set D to a new value.
  D = 2.*B + 3.*C +2.*D;
  std::cout << D << std::endl;
  // get a pointer to the begining of the internal 1D array that store the value of D. (this is usefull to pass the raw data to a C function or fortran function or to work directly on the raw data, without using the (i,j) operator It is often used in the dmatrix_denseCM member function)
  double *a = D.data();
  // print the content of the internal array. note that the componnent of D are internally stored Column by Column. This is Called the Columns Major Format.
  int m = D.getNbLines();
  int n = D.getNbColumns();
  for(int k = 0; k < m*n ; ++k ){
    std::cout << a[k] << " ";
  }
  std::cout << std::endl;
  // the term i,j of the matrix is therefore at position k=i+j*m is the data:
  int i = 1, j=1;
  if((i >= m) || (j>=n)) std::cout << "matrix D too small to continue examples. stop."<< std::endl; return 1;
  int k = i+j*m;
  std::cout <<"D(" << i << "," << j <<") = " <<  D(i,j) << std::endl;
  std::cout << "a[" << k << "] = " << a[k] << std::endl;
  // the k th term in the data correspond to the following i,j term :
  k = 3;
  i = k%m;
  j = k/m;
  std::cout <<"D(" << i << "," << j <<") = " <<  D(i,j) << std::endl;
  std::cout << "a[" << k << "] = " << a[k] << std::endl;
  // iterate over the entry of D, column by column
  for( auto v: D) std::cout << v << " ";
  std::cout << std::endl;

  // try to open the file "data.mat", that contain an example of a matrix stored in file, readable by the dmatrix_denseCM class.
  std::ifstream inputfile("data.mat");
  // check if the file was found.
  if(!inputfile.good()){
     std::cout << "Can't open file data.mat" << std::endl;
     throw;
  }
  // try to read the input file and store the result in matrix D 
  inputfile >> D;
  // print the new D on screen
  std::cout << "D as read from file data.mat " << std::endl;  
  std::cout << D << std::endl;


  // Generating random numbers :
  // First set a random number generator
  std::random_device rd;
  std::mt19937 gen(rd());  // this is the generator
  std::uniform_real_distribution<> dis(0., 1.); // this is to have an uniform distribution between 0. and 1.
  // get and print 10 random number between 0. and 1., using our generator
  for (int i = 0; i< 10; ++i){
    std::cout << dis(gen) << std::endl;
  }
  // Fill up the Matrix D with random number
  for (auto &v : D) v = dis(gen);
  std::cout << D << std::endl;
}
