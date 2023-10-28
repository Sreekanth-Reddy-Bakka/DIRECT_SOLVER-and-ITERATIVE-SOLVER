#include <math.h>
#include "blas_lapack_def.h"
#include "dLU_denseCM.h"

int factorBasic(int n,  double *a, int LDA ){
  std::cout <<" Error factorBasic Not implemented " << __FILE__ << __LINE__ << std::endl; 
  throw;
}

int factorL2(int n,  double *a, int LDA ){
   std::cout <<" Error factorL2 Not implemented " << __FILE__ << __LINE__ << std::endl; 
  throw;
}

int factorL3(int r, int n,  double *a, int LDA ){
  std::cout <<" Error factorL2 Not implemented " << __FILE__ << __LINE__ << std::endl; 
  throw;
}
 

void factorlapack( int m, double *a , int LDA, int *ipiv){
  int info;
  dgetrf_(m,m, a, LDA, ipiv, info);
  
  if (info != 0){
    std::cout << "Factorisation using lapack dgetrf2_ failed ! " << std::endl;
    std::cout << "dgetrf2 info " << info << std::endl;
    throw;
  }
}

dLU_denseCM::dLU_denseCM(const dmatrix_denseCM &A, factorpolicy fpol ):n(A.getNbLines()), LDA(n), a(new double[n*n]), ipiv(new int[n]) {
  assert(A.getNbLines()==A.getNbColumns());
  std::copy(A.begin(), A.end(), a );
  switch (fpol.impl){
  case (factorpolicy::Basic):{
    factorBasic(n, a, LDA);
    for (int i = 0; i < n ;++i) ipiv[i] = i+1;
    break;
  }
  case(factorpolicy::L2):{
    factorL2(n, a, LDA);
    for (int i = 0; i < n ;++i) ipiv[i] = i+1;
    break;
  }
  case(factorpolicy::L3):{
    int blocksize = fpol.block_size;
    factorL3(blocksize, n, a, LDA);
    for (int i = 0; i < n ;++i) ipiv[i] = i+1;
    break;
  }
  case(factorpolicy::Lapack):{
    factorlapack(n,a,LDA,ipiv);
    break;
  }
  }
}

dmatrix_denseCM dLU_denseCM::getL() const{
  dmatrix_denseCM L(n,n,0.);   
  for (int j  = 0; j < n; ++j){
    L(j,j) = 1.;
    for (int i=j+1; i < n; ++i){
	L(i,j) = a[i+j*LDA];
    }
  }
  return L;
}

dmatrix_denseCM dLU_denseCM::getU() const{
  dmatrix_denseCM U(n,n,0.);   
  for (int j  = 0; j < n; ++j){
    for (int i=0; i < j+1; ++i){
      U(i,j) = a[i+j*LDA];
    }
  }
  return U;
}
  
dmatrix_denseCM dLU_denseCM::getP() const{
  dmatrix_denseCM P(n,n,0.);
  for (int i =0; i < n; ++i){
    P(i,i) = 1;
  }
  for (int i =0; i < n; ++i){
    int j = ipiv[i] -1;
    if(i!=j){
      for (int k = 0; k < n; ++k) std::swap(P(i,k), P(j,k));
    }
  }
  return P;
}




dLU_denseCM::~dLU_denseCM(){
  delete[] a;
  delete[] ipiv;
}


dmatrix_denseCM dLU_denseCM::solve(const dmatrix_denseCM &B) const {
    int nrhs = B.getNbColumns();
    int nequ = B.getNbLines();
    assert(nequ == n);
    dmatrix_denseCM X(B);
    for (int i = 0; i < nequ; ++i){
      int j = ipiv[i]-1; 
      if(i!=j){ 
	dswap_(nrhs, X.data()+i, nequ, X.data()+j, nequ );
      }
    }
    dtrsm_('L', 'L', 'N', 'U', n, nrhs, 1., a, n, X.data(), nequ );
    dtrsm_('L', 'U', 'N', 'N', n, nrhs, 1., a, n,  X.data(), nequ );
    return X; 
  }


 
