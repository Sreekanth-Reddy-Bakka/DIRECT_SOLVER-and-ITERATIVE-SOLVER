#ifndef _dLLT_denseCM_
#define _dLLT_denseCM_
#include "dmatrix_denseCM.h"
class dLLT_denseCM{
 public :
  // 
  dLLT_denseCM(const dmatrix_denseCM &A );
  ~dLLT_denseCM();
  /// This solve the linear Ax =B, using the factorisation. The returned matrix is X.

  dmatrix_denseCM solve(const dmatrix_denseCM &B) const;
  dmatrix_denseCM getL() const;
  /// This function return a dmatrix_denseCM containing L ... This is mostly for illustration purpose.
  private:
  int n;
  int LDA;
  double *a;
 };

#endif
