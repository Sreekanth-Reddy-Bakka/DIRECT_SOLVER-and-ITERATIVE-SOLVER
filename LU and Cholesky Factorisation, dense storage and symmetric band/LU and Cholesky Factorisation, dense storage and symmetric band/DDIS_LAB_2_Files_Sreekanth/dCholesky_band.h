#ifndef _dCholeskyBand_
#define _dCholesky_band_
#include "dsquarematrix_symband.h"
#include "dmatrix_denseCM.h"


/// This class is meant to store a Cholesky (LLT ) factorisation of a symmetric ositive matrix in band format.
class dCholesky_band{
 public:
  /// Construction do the factorisation of the input dsquarematrixinsymetric band format, and store the result in the object
  dCholesky_band(const dsquarematrix_symband &in );
  ~dCholesky_band();
  /// Solve using the factorisation for the right hand side given as a densematrix in column major format.
  dmatrix_denseCM solve(const dmatrix_denseCM &B) const;
 private :
  int m;
  int lb;
  double *a;
};
#endif
