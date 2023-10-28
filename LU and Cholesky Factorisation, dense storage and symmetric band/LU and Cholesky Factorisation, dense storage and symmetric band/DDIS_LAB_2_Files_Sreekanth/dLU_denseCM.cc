#include <math.h>
#include "blas_lapack_def.h"
#include "dLU_denseCM.h"

int factorBasic(int n,  double *a, int LDA )
{
	int k = 0;
	double piv = a[0];
	const double min_piv= 1e-6;
	
	while((std::fabs(piv) > min_piv) && (k < n-1)) 
	{
		for (int i = k+1; i < n; ++i)
		{
			a[i + k*LDA] /= piv;
		}
		
		for (int i = k+1; i < n; ++i)
		{
			for (int j = k+1; j < n; ++j)
			{
				a[i + j*LDA] -= a[i + k*LDA] * a[k + j*LDA];
			}
		}
		k += 1;
		piv = a[k + k*LDA];
	}
	
	if (std::fabs(piv ) <= min_piv)
	{
		std::cout << "NULL PIVOT IN dLU1" << piv << std::endl;
		return 0;
	}
	return 1;
}

int factorL2(int n,  double *a, int LDA )
{
	int k = 0;
	double piv = a[0];
	const double min_piv= 1e-6;
	
	while((std::fabs(piv) > min_piv) && (k < n-1))
	{
		dscal_(n-(k+1), 1./piv, a+(k+1+k*LDA), 1);
		dger_(n - (k+1), n - (k+1), -1., a+(k+1 + k*LDA), 1, a+(k+(k+1)*LDA), LDA, a+(k+1+(k+1)*LDA), LDA);
		k += 1;
		piv = a[k+k*LDA];
	}
	if (std::fabs(piv ) <= min_piv)
	{
		std::cout << "NULL PIVOT IN dLU2" << piv << __FILE__ << ":" << __LINE__ << std::endl;
		return 0;
	}
	return 1;
}

// r is the block size:
int factorL3(int r, int n,  double *a, int LDA )
{
	int  l = 0;
	while (l < n)
	{
		int m = std::min(n, l+r); // m = 3
		int bsize = m - l; // bsize = 3
		int success = factorL2( bsize, a+(l + l*LDA), LDA);
		if(!success)
		{
			std::cout << "CAN'T FACTORIZE ONE BLOCK" << std::endl;
			return 0;
		}
		
		dtrsm_('L', 'L', 'N', 'U', bsize, n-m, 1., a+(l+l*LDA), LDA, a+l+m*LDA, LDA);
		dtrsm_('R', 'U', 'N', 'N', n-m, bsize, 1., a+(l+l*LDA), LDA, a+m+l*LDA, LDA);
		
		dgemm_('N', 'N', n-m, n-m, bsize, -1, a + m + l*LDA, LDA, a + l + m*LDA, LDA, 1., a+m+m*LDA, LDA);
		l = m;		
	}
	return 1;
}

void factorlapack( int m, double *a , int LDA, int *ipiv)
{
	int info;
	dgetrf2_(m,m, a, LDA, ipiv, info);
	if (info != 0)
	{
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
