#include <algorithm>
#include <assert.h>
#include <iostream>
#include <functional>
#include <limits>
#include <math.h>
#include<numeric>
#include "dmatrix_denseCM.h"
#include "mmio.h"

#include "blas_lapack_def.h"

/**************************************!
   YOU HAVE TO IMPLEMENT THIS FUNCTION !!!
   it must return the upper bandwidth of A
 **************************************/  
int computeBandwidthUp(const dmatrix_denseCM &A )
{
	int m = A.getNbLines();
	int n = A.getNbColumns();
	int upperBandwidth = 0;
	for (int i = 0; i < m; i++) 
	{
        	for (int j = 0; j < n; j++) 
        	{
        		if (i < j)
        		{
				if (A(i,j) != 0) 
				{
					int distance = j - i;				
					if (distance > upperBandwidth) 
					{
						upperBandwidth = distance; 
					}
				}
			}
        	}
        }
        return upperBandwidth;
};

/**************************************!
   YOU HAVE TO IMPLEMENT THIS FUNCTION !!!
   it must return the upper bandwidth of A
 **************************************/  
int  computeBandwidthDown(const dmatrix_denseCM &A )
{
	int m = A.getNbLines();
	int n = A.getNbColumns();
	int lowerBandwidth = 0;
	for (int i = 0; i < m; i++) 
	{
        	for (int j = 0; j < n; j++) 
        	{
			if (i > j)
        		{
				if (A(i,j) != 0) 
				{
					int distance = i - j;				
					if (distance > lowerBandwidth) 
					{
						lowerBandwidth = distance; 
		        		}
		        	}
		        }
        	}
        }
        return lowerBandwidth;
};

dmatrix_denseCM::dmatrix_denseCM(size_t _m, size_t _n, double val):m(_m), n(_n), a(new double [m*n]) 
{
	std::fill(a, a+m*n, val);
}

dmatrix_denseCM::dmatrix_denseCM(size_t _m, size_t _n):m(_m), n(_n), a(new double [m*n]) 
{

}

dmatrix_denseCM::dmatrix_denseCM(const dmatrix_denseCM &in):m(in.m), n(in.n), a(new double [m*n]) 
{
	std::copy(in.begin(), in.end(), a);
}

dmatrix_denseCM::dmatrix_denseCM(size_t _m, size_t _n, const double *a, storage s):m(_m), n(_n), a(new double[m*n])
{
	setData(a, a+m*n, s);
}

dmatrix_denseCM& dmatrix_denseCM::operator = (const dmatrix_denseCM &in)
{
	reallocate(in.m, in.n);
	std::copy(in.begin(), in.end(), data());
	return *this;
}

dmatrix_denseCM::~dmatrix_denseCM()
{
	delete []a;
}

int dmatrix_denseCM::getNbLines() const 
{
	return m;
}

int dmatrix_denseCM::getNbColumns() const 
{
	return n;
}

const double * dmatrix_denseCM::data() const 
{
	return a;
}

double * dmatrix_denseCM::data()
{
	return a;
}

void dmatrix_denseCM::setData(const double *begin, const double *end, storage s )
{
	assert (std::distance(begin, end) <= m*n);
	switch (s)
	{
		case CM:
		{
			std::copy(begin, end, data());
			return;
		}
		case RM:
		{
			int ndata = std::distance(begin, end);
			int nfull_row = ndata/n;
			for (int i = 0; i < nfull_row; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					(*this)(i,j) = *(begin + i*n +j);
				}
			}
			for (int j = 0; j < ndata%n; ++j)
			{
				(*this)(nfull_row,j) = *(begin + nfull_row*n +j);
			}
			return;
		}
	}
}

double *dmatrix_denseCM::begin()
{
	return a;
}

double *dmatrix_denseCM::end()
{
	return a+m*n;
} 
 
const double *dmatrix_denseCM::begin() const
{
	return a;
}

const double *dmatrix_denseCM::end()const 
{
	return a+m*n;
} 

double &dmatrix_denseCM::operator()(size_t i, size_t j)
{
	assert(i<m && j <n );
	return *(a+i+m*j);
}

const double &dmatrix_denseCM::operator()(size_t i, size_t j) const
{
	assert(i<m && j <n );
	return *(a+i+m*j);
}


dmatrix_denseCM &dmatrix_denseCM::operator+=(const dmatrix_denseCM &r)
{
	assert ((r.m*r.n) == (m*n));
	std::transform(begin(), end(), r.begin(), begin(), std::plus<double>() );
	return *this;
}
  
dmatrix_denseCM &dmatrix_denseCM::operator-=(const dmatrix_denseCM &r)
{
	assert ((r.m*r.n) == (m*n));
	std::transform(begin(), end(), r.begin(), begin(), std::minus<double>() );
	return *this;
}

dmatrix_denseCM &dmatrix_denseCM::operator*=(double a)
{
	std::for_each(begin(), end(), [&a](double &in ){in*=a;});
	return *this;
}

dmatrix_denseCM &dmatrix_denseCM::operator/=(double a)
{
	std::for_each(begin(), end(), [&a](double &in ){in/=a;});
	return *this;
}

/// return sqrt (sum a_ij^2)
double dmatrix_denseCM::norm2() const{
  struct op{
    double operator()(double sum, double a){
      return sum + a*a;
    }
  };
  
  double sum = std::accumulate(begin(), end(), 0., op() );
  return sqrt(sum);
} 

void dmatrix_denseCM::reallocate(int _m, int _n){
  if(m*n != _m*_n) {
    delete [] a;
    a = nullptr;
  }
  m = _m;
  n = _n;
  if(!a) a = new double[m*n];
  return;
}

dmatrix_denseCM operator-(const dmatrix_denseCM &l, const dmatrix_denseCM &r ){
  dmatrix_denseCM res(l);
  res-=r;
  return res; 
}

dmatrix_denseCM operator+(const dmatrix_denseCM &l, const dmatrix_denseCM &r ){
  dmatrix_denseCM res(l);
  res+=r;
  return res; 
}
dmatrix_denseCM operator*(const dmatrix_denseCM &l, double a ){
  dmatrix_denseCM res(l);
  res*=a;
  return res; 
}

dmatrix_denseCM operator*(double a, const dmatrix_denseCM &r){
  return r*a;
}



dmatrix_denseCM mulijk(const dmatrix_denseCM &A, const dmatrix_denseCM &B){
  assert(A.getNbColumns() == B .getNbLines());
  const size_t m = A.getNbLines();
  const size_t n = B.getNbColumns();
  const size_t o = B.getNbLines();
  dmatrix_denseCM C(m,n,0.);
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      for (size_t k = 0; k < o; ++k)
	C(i,j) += A(i,k)*B(k,j);
  return C;
}

dmatrix_denseCM muljki(const dmatrix_denseCM &A, const dmatrix_denseCM &B){
  assert(A.getNbColumns() == B .getNbLines());
  const size_t m = A.getNbLines();
  const size_t n = B.getNbColumns();
  const size_t o = B.getNbLines();
  dmatrix_denseCM C(m,n,0.);
  
  for (size_t j = 0; j < n; ++j)
    for (size_t k = 0; k < o; ++k){
      //      const double Bkj = B(k,j);
      //#pragma omp parallel for numthreads(2)
      for (size_t i = 0; i < m; ++i)
	C(i,j) += A(i,k)*B(k,j);
    }
  return C;
}


dmatrix_denseCM muldgemm(const dmatrix_denseCM &A, const dmatrix_denseCM &B ){
  assert(A.getNbColumns() == B .getNbLines());
  const int  m = A.getNbLines();
  const int  n = B.getNbColumns();
  const int  k = B.getNbLines();
  dmatrix_denseCM C(m,n);
  dgemm_('N', 'N', m, n, k, 1., A.data(), m, B.data(), k, 0., C.data(), m );
  return C;
}

double norm2(const dmatrix_denseCM &A )
{
	return A.norm2();
}
 
std::ostream& operator<<(std::ostream& out, const dmatrix_denseCM& A)
{
  int m = A.getNbLines();
  int n = A.getNbColumns();
  out << m <<" " << n  << std::endl;
  for (size_t i =0; i < m; ++i) {
    for (size_t j =0; j < n; ++j) {
      out <<  A(i,j) << " ";
    }
    out << std::endl;
  }
  return out;
}

std::istream& operator>>(std::istream& in, dmatrix_denseCM& A)
{
  bool matrixmarket = false;
  if (&in == &std::cin){
    std::cout << "input number of line, number of column, then all the terms line by line"<< std::endl;
  }
  char next;
  next =in.peek();
  if (next == '%') {
    in.get();
    if (in.peek() == '%'){
      //We are working with a matrix market file !
      matrixmarket = true;
    }
    in.unget();
  }
  if (matrixmarket) {
    readMatrixMarket(in,A);
    return in;
  }
  readMyFormat(in, A);
  return in;
}

void readMyFormat(std::istream &in, dmatrix_denseCM &A){
  char next;
  next = in.peek();
  while(next=='#'||next=='%'){
    in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    next = in.peek();
  }
  int _m, _n;
  in >> _m;
  if(!in.good()) {std::cout << "I could not read an Integer for the number of lines" << std::endl; throw;}
  in >> _n;
  if(!in.good()) {std::cout << "I could not read an Integer for the number of columns" << std::endl; throw;}
  double *ain = new double[_m*_n];
  for (int i = 0; i < _m*_n; ++i) {
    in >> ain[i];
    if(in.fail()) {
      int line   = i/_m;
      int column = i%_m;
      std::cout << "I could not read the " << i << "th entry of the matrix (" << line<< "," << column << ") " << std::endl; throw;
    }
  }
  A.reallocate(_m,_n);
  A.setData(ain, ain+_m*_n, RM);
  
  delete []ain;
}


void readMatrixMarket(std::istream &in, dmatrix_denseCM &A){
  int ret_code;
  MM_typecode matcode;
  int M, N, nz; 
  if (mm_read_banner(in, &matcode) != 0)
    {
      std::cout << "Could not process Matrix Market banner" <<  std::endl;
      throw;
    }
  
  if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
      mm_is_sparse(matcode)  )
    {
      std::cout << "Sorry, this application does not support ";
      std::cout << "Market Market type: ["  << mm_typecode_to_str(matcode);
      std::cout <<   "]" << std::endl;
      exit(1);
    }
  if mm_is_sparse(matcode){
      int M, N, nz;
      if ((ret_code = mm_read_mtx_crd_size(in, &M, &N, &nz)) !=0){
	printf("ARG");
	exit(1);
      }
      A.reallocate(M, N);
      std::fill(A.data(), A.data()+M*N, 0.);
      int i,j;
      double val;
      for (int k=0; k<nz; k++)
	{
	   in >> i >> j >> val;
	   A(i-1,j-1) = val; 
	   if (mm_is_symmetric( matcode)) A(j-1,i-1) = val; 
	}
    }
}



dmatrix_denseCM transpose(const dmatrix_denseCM &A){
  int m =  A.getNbColumns();
  int n  = A.getNbLines(); 
  dmatrix_denseCM res(m,n);
  for (int i =0 ; i <m ; ++i) {
    for (int j = 0; j < n; ++j)
      res(i,j) = A(j,i);
  }
  return res;
}

dmatrix_denseCM operator*(const dmatrix_denseCM &A, const dmatrix_denseCM &B)
{
	return muldgemm(A,B);	
}

void permutlines(dmatrix_denseCM &A, int *ipiv)
{
	int m = A.getNbLines();
	int n = A.getNbColumns();
	for (int i = m-1; i >= 0 ; --i)
	{
		int j = ipiv[i] - 1;
		if (i != j)
		{
			for (int k = 0; k < n; ++k)
			std::swap(A(i,k), A(j,k));
		}
	}
}
