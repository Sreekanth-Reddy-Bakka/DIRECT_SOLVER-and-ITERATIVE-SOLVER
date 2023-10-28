#ifndef _dLLT_CCS_
#define _dLLT_CCS_
#include "dmatrix_CCS.h"
#include "graph.h"

/// This class Stores a cholesky factorisation, and can solve an associated linear system.
class dLLT_CCS{
 public:
  /// Cholesky factorisation constructor. A must be Symetric positive definite
  /*! A = LLT =UTU */
  dLLT_CCS( const dmatrix_CCS &A);
  /// return LT in CSS format
  dmatrix_CCS getLT() const;
  /// return L  in CSS format
  dmatrix_CCS getL()  const;
  /// return the solution X of AX = B. 
  /*!  B must be a dmatrix_denseCM, must have same number of lines that A
       Each column of B represent one rightand side to solve for 
       The function return the matrix X, solution of the linear systems*/
  dmatrix_denseCM solve( const dmatrix_denseCM &B);
 private :
  dmatrix_CCS LT;
  int n;
  graph LpLTgraph (const graph & gA);
  graph colgraphL (const graph & gLpLT);
  graph colgraphLT (const graph & gLpLT);

};
#endif
