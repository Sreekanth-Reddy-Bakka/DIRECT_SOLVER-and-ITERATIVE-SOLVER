#include "dmatrix_denseCM.h"

class dmatrix_view_denseCM{
 dmatrix_view_denseCM( dmatrix_denseCM &_A, size_t _i, size_t _j, size_t _m, size_t _n ):A(_A), si(_i), sj(_j), m(_m), n(_n){
    
  }
  double & operator()(size_t _i, size_t _j){
    return A(si+m, sj+n);
  }
  double *getData(){
    return &A(si,sj);
  }
  const double *getData() const{
    return &A(si,sj);
  }
  size_t LD() const {
    return A.LD;
  }
  dmatrix_denseCM &A;
  size_t si, sj, m, n;
};

void dgemm( const double &alpha, const dmatrix_view_denseCM &A, const dmatrix_view_denseCM &B, const double &beta, dmatrix_view_denseCM &C){
  assert(A.n == B.m);
  dgemm_( 'N', 'N', A.m, B.n, A.n, alpha, A.data(), A.LD, B.data(), B.LD, beta, C.data(), C.LD);
}


