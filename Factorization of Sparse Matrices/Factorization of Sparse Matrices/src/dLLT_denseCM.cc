#include <math.h>
#include "blas_lapack_def.h"
#include "dLLT_denseCM.h"


dLLT_denseCM::dLLT_denseCM(const dmatrix_denseCM &A ):n(A.getNbLines()), LDA(n), a(new double[n*n]) {
  // good version per lines
  if(1){
    std::fill(a,a+n*n,0.);
    for (int j = 0; j < n; ++j){
      double tmp  = A(j,j);   
      for (int k = 0; k < j; ++k) tmp -= a[j+k*LDA]*a[j+k*LDA];
      a[j+j*LDA] = sqrt(tmp);
      for (int i = j+1 ; i < n; ++i){
	tmp = A(i,j);
	for (int k = 0; k < j; ++k) tmp -= a[i+k*LDA]*a[j+k*LDA];
	a[i+j*LDA] = tmp/a[j+j*LDA];
      }
    }
  }
  //column version compute u ...
  if(0){
    std::vector<double> u(n*n);
    std::fill(u.begin(),u.end(),0.);
    for (int j = 0; j < n; ++j){
      double tmp  = A(j,j);   
      for (int k = 0; k < j; ++k) tmp -= u[k+j*LDA]*u[k+j*LDA];
      u[j+j*LDA] = sqrt(tmp);
      for (int i = j+1 ; i < n; ++i){
	tmp = A(i,j);
	for (int k = 0; k < j; ++k) tmp -= u[k+i*LDA]*u[k+j*LDA];
	u[j+i*LDA] = tmp/u[j+j*LDA];
      }
    }
    
    // then transpose the result to get L
    std::fill(a,a+n*n,0.);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j <= i; ++j)
	a[i+j*LDA] = u[j+i*LDA];
  }
}

dLLT_denseCM::~dLLT_denseCM(){delete []a;}
  /// This solve the linear Ax =B, using the factorisation. The returned matrix is X.

dmatrix_denseCM dLLT_denseCM::getL() const{
  dmatrix_denseCM L(n,n);
  std::copy(a,a+n*n, L.data());
  return L;
}

dmatrix_denseCM dLLT_denseCM::solve(const dmatrix_denseCM &B) const{
  int nrhs = B.getNbColumns();
  int nequ = B.getNbLines();
  assert(nequ == n);
  dmatrix_denseCM X(B);
  dtrsm_('L', 'L', 'N', 'N', n, nrhs, 1., a, n, X.data(), nequ );
  dtrsm_('L', 'L', 'T', 'N', n, nrhs, 1., a, n,  X.data(), nequ );
    return X; 
    
}
