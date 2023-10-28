#include "dLLT_CCS.h"
#include "dmatrix_CCS.h"
#include "graph.h"
#include <iomanip>
#include <math.h>

dLLT_CCS::dLLT_CCS( const dmatrix_CCS &A):n(A.m){
    assert( n == A.n);
    graph gA = buildSymGraph(A);
    graph gLT = colgraphLT(LpLTgraph(gA));
    LT = dmatrix_CCS (gLT, 0.); 
    for (int j = 0; j < n; ++j){
      double tmp  = A(j,j);
      // look at col j computer inner product of colmnj 
      for (int k = LT.columnptr[j]; k < LT.columnptr[j+1]; ++k) tmp -= LT.a[k]*LT.a[k];
      LT(j,j) = sqrt(tmp);
      // update lines j col i
      for (int i = j+1 ; i < n; ++i){
	int kk = LT.columnptr[i];
	for ( ; kk <LT.columnptr[i+1]; ++kk){
	  if (j == LT.lineindex[kk]) break;
	}
	if (kk < LT.columnptr[i+1]){
	  tmp = A(i,j);
	  // assuming both line are sorted ... which is the case in present context (graphs are sorted !)
	  int ki = LT.columnptr[i];
	  int eki = LT.columnptr[i+1];
	  int rki = LT.lineindex[ki];	  
	  for (int kj = LT.columnptr[j]; kj < LT.columnptr[j+1] , ki < eki  ; ++kj) {
	    int rkj = LT.lineindex[kj];
	    for( ;  ki < eki; ++ki) {
	      rki = LT.lineindex[ki];
	      if (rki >= rkj) break;
	    }
	    if (rki == rkj){
	      tmp -= LT.a[ki]*LT.a[kj];
	    }
	  }
	  LT(j,i) = tmp/LT(j,j);
	}
      }
    }  
  }



graph dLLT_CCS::LpLTgraph (const graph & gA){
  int n = gA.nbNodes();
  graph g(gA);
  for (int k = 0; k < n; ++k){
    auto N = g.getNeighbour(k);
    for (int i : N){
      if ( i > k){
	for (int j : N){
	  if (j > k){
	    g.addEdge(i,j);
	  }
	}
      }
    }
  }
  return g;
}

graph dLLT_CCS::colgraphL (const graph & gLpLT){
  graph gL;
  gL.addNodes(n);
  for (int j = 0; j < n; ++j){
    for ( int i : gLpLT.getNeighbour(j )){
      if (i >= j) gL.addEdge(j,i);
    }
  }
  return gL;
}

graph dLLT_CCS::colgraphLT (const graph & gLpLT){
  graph gLT;
  gLT.addNodes(n);
  for (int j = 0; j < n; ++j){
    for ( int i : gLpLT.getNeighbour(j )){
      if (i <= j) gLT.addEdge(j,i);
    }
  }
  return gLT;
}
  
dmatrix_CCS dLLT_CCS::getLT()const {
  return LT;
}

dmatrix_CCS dLLT_CCS::getL()const {
  return transpose(LT);
}

dmatrix_denseCM dLLT_CCS::solve( const dmatrix_denseCM &B){
  dmatrix_denseCM X(B);
  return solveU(LT, solveUT(LT,X));
}


