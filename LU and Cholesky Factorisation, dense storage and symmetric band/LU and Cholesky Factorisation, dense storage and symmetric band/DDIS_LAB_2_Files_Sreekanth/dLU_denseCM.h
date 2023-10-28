#ifndef _dLU_denseCM_
#define _dLU_denseCM_
#include "dmatrix_denseCM.h"

/// This function is meant to factorise a square matrix in Column major storage
/* This function must implement the basic algorithm: straight loops. you have to implement it.
   on entry :
   a point to the begining of the matrix storage or to the submatrix you want to factorise.
   n is the number of line (and column) of your block.
   LDA, is the number of term you need to jump to get to the beggining of the next column of your submatrix (if your submatrix is the full matrix LDA =n).
*/
int factorBasic(int n,  double *a, int LDA );
/// This function is meant to factorise a square matrix in Column major storage
/* This function must implement the Level 2 blas version algorithm: inner loops are replaced by call to level1blas function (dger and dscal)
   on entry :
   a point to the begining of the matrix storage or to the submatrix you want to factorise.
   n is the number of line (and column) of your block.
   LDA, is the number of term you need to jump to get to the beggining of the next column of your submatrix (if your submatrix is the full matrix LDA =n).
*/
int factorL2(int n,  double *a, int LDA );
/// This function is meant to factorise a matrix in Column major storage
/* This function must implement the Level 3 blas version algorithm: the block factoring version
   on entry :
   a point to the begining of the matrix storage or to the submatrix you want to factorise.
   n is the number of line (and column) of your block.
   LDA, is the number of term you need to jump to get to the beggining of the next column of your submatrix (if your submatrix is the full matrix LDA =n).   
   bs is the block size.
*/
int factorL3(int bs, int n,  double *a, int LDA );

/// This function do a factorisation with pivoting, using Lapack. It is already implemented and you don't have to do any thing here.
void factorlapack( int m, double *a , int LDA, int *ipiv);

class factorpolicy{
 public :
  enum factorimpl{Basic, L2, L3, Lapack};
 factorpolicy(factorimpl _impl=Basic, int  _block_size=0):impl(_impl), block_size(_block_size){} ;
  const factorimpl impl;
  const int block_size;
};

/// This class compute and store an LU factorisation of a dense Matrix.
/// The factors are computed at construction, using one of the implemeted algorithm, depending on the parameters passed at construction.
class dLU_denseCM{
 public :
  // 
  dLU_denseCM(const dmatrix_denseCM &A, factorpolicy fpol );
  ~dLU_denseCM();
  /// This solve the linear Ax =B, using the factorisation. The returned matrix is X.
  dmatrix_denseCM solve(const dmatrix_denseCM &B) const;
  /// This function return a dmatrix_denseCM containing L ... This is mostly for illustration purpose.
  /// No need to use it in serious code.
  dmatrix_denseCM getL() const;
  /// This function return a dmatrix_denseCM containing U ... This is mostly for illustration purpose.
  /// No need to use it in serious code.
  dmatrix_denseCM getU() const;
  /// This function return a dmatrix_denseCM containing P ... This is mostly for illustration purpose.
  /// No need to use it in serious code.
  dmatrix_denseCM getP() const;
 private:
  int n;
  int LDA;
  double *a;
  int *ipiv;
};

#endif
