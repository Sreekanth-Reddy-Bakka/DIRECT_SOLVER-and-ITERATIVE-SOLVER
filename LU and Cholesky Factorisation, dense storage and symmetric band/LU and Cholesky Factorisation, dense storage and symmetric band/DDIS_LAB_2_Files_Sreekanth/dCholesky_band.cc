#include "dCholesky_band.h"
#include "blas_lapack_def.h"

dCholesky_band::dCholesky_band(const dsquarematrix_symband &in ):m(in.getNbLines()), lb(in.getBandwidth()), a(new double[m*(lb+1)]){
  std::copy(in.begin(), in.end(), a);
  int info;
  dpbtrf_('U', m, lb, a, lb+1, info );
  if(info != 0){
    std::cout << "Cholesky Factorisation failed ! Some minor are not postive definite " << __FILE__ << " " << __LINE__ <<std::endl;
    throw;
  }
}

dCholesky_band::~dCholesky_band(){delete []a;}

dmatrix_denseCM dCholesky_band::solve(const dmatrix_denseCM &B) const{
  int nrhs = B.getNbColumns();
  int nequ = B.getNbLines();
  assert(nequ == m);
  dmatrix_denseCM X(B);
  int info;
  dpbtrs_('U',m,lb, nrhs, a, lb+1, X.data(), nequ, info);
  if (info != 0){
    std::cout << "Error calling lapack dpbtrs, argument " << -info << " has an illegal value " << __FILE__<< " "<<__LINE__<< std::endl;
    throw;
  }
  return X;
}

