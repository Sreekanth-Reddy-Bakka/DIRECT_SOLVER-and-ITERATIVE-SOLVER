#include <algorithm>
#include <assert.h>
#include <iostream>
#include <functional>
#include <limits>
#include <math.h>
#include<numeric>
#include "dmatrix_denseCM.h"
#include "mmio.h"

extern "C"{
  void dgemm_( const char &TRANSA, const char &TRANSB, const int &M, const int  &N, const int & K,  const double & ALPHA,  const double *A, const int & LDA, const double *B,  const int &LDB, const double &BETA, double *C, const int & LDC);
}


dmatrix_denseCM::dmatrix_denseCM(size_t _m, size_t _n, double val):m(_m), n(_n), a(new double [m*n]) {
  std::fill(a, a+m*n, val);
}

dmatrix_denseCM::dmatrix_denseCM(size_t _m, size_t _n):m(_m), n(_n), a(new double [m*n]) {
}

dmatrix_denseCM::dmatrix_denseCM(const dmatrix_denseCM &in):m(in.m), n(in.n), a(new double [m*n]) {
  std::copy(in.begin(), in.end(), a);
}

dmatrix_denseCM::dmatrix_denseCM(size_t _m, size_t _n, const double *a, storage s):m(_m), n(_n), a(new double[m*n]){
  setData(a, a+m*n, s);
}

dmatrix_denseCM& dmatrix_denseCM::operator = (const dmatrix_denseCM &in){
  reallocate(in.m, in.n);
  std::copy(in.begin(), in.end(), data());
  return *this;
}

dmatrix_denseCM::~dmatrix_denseCM(){
  delete []a;
}

int dmatrix_denseCM::getNbLines() const {
  return m;
}

int dmatrix_denseCM::getNbColumns() const {  
  return n;
}

const double * dmatrix_denseCM::data() const {
  return a;
}

double * dmatrix_denseCM::data(){
  return a;
}

void dmatrix_denseCM::setData(const double *begin, const double *end, storage s ){
    assert (std::distance(begin, end) <= m*n);
    switch (s){
    case CM:{
      std::copy(begin, end, data());
      return;
    }
    case RM:{
      int ndata = std::distance(begin, end);
      int nfull_row = ndata/n;
      for (int i = 0; i < nfull_row; ++i){
	for (int j = 0; j < n; ++j){
	  (*this)(i,j) = *(begin + i*n +j);
	}
      }
      for (int j = 0; j < ndata%n; ++j){
	(*this)(nfull_row,j) = *(begin + nfull_row*n +j);
      }
      return;
    }
    }
}

double *dmatrix_denseCM::begin(){
  return a;
}

double *dmatrix_denseCM::end(){
  return a+m*n;
} 
 
const double *dmatrix_denseCM::begin() const{
  return a;
}

const double *dmatrix_denseCM::end()const {
  return a+m*n;
} 

double &dmatrix_denseCM::operator()(size_t i, size_t j){
  assert(i<m && j <n );
  return *(a+i+m*j);
}

const double &dmatrix_denseCM::operator()(size_t i, size_t j) const{
  assert(i<m && j <n );
  return *(a+i+m*j);
}


dmatrix_denseCM &dmatrix_denseCM::operator+=(const dmatrix_denseCM &r){
  assert ((r.m*r.n) == (m*n));
  std::transform(begin(), end(), r.begin(), begin(), std::plus<double>() );
  return *this;
}
  
dmatrix_denseCM &dmatrix_denseCM::operator-=(const dmatrix_denseCM &r){
  assert ((r.m*r.n) == (m*n));
  std::transform(begin(), end(), r.begin(), begin(), std::minus<double>() );
  return *this;
}

dmatrix_denseCM &dmatrix_denseCM::operator*=(double a){
  std::for_each(begin(), end(), [&a](double &in ){in*=a;});
  return *this;
}

dmatrix_denseCM &dmatrix_denseCM::operator/=(double a){
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
  
  double sum = std::accumulate(begin(), end(), 0, op() );
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

bool operator==(const dmatrix_denseCM &lhs, const dmatrix_denseCM &rhs)
{
	const size_t M = lhs.getNbLines();
	const size_t K = lhs.getNbColumns();
	const size_t N = lhs.getNbColumns();	
	
	for (size_t i = 0; i < M; ++i) 
	{
		for (size_t j = 0; j < N; ++j) 
		{
			if ((lhs(i, j) - rhs(i, j)) > 0.00000000000000001)
			{
				return false; 
			}
		}
	}

	return true;
	
}

dmatrix_denseCM mulV1(const dmatrix_denseCM &A, const dmatrix_denseCM &B){
  // get the size data from A and B
  const size_t M = A.getNbLines();
  const size_t K = A.getNbColumns();
  const size_t N = B.getNbColumns();
  if (K != B.getNbLines()){
    // The code will stop if there is a missmatch on the size of A and B
      std::cout <<"error can't compute product of 2 matrices : dimension missmatch." << std::endl;
      throw;
  }
  // The size are correct, the matrix C that will store the result is created
  dmatrix_denseCM C(M,N,0.);
  
  const double * a = A.data();
  const double * b = B.data();
  double * c = C.data();
  
  // This is where you need to work : implement your matrix multiplication below, before return C

  for(int j = 0; j < N; ++j)
  {
  	for(int k = 0; k < K; ++k)
  	{
  		for(int i = 0; i < M; ++i)
  		{
  			*(c+i+j*M) += (*(a+i+k*M)) * (*(b+k+j*K));
  		}
  	}
  }

  //std::cout << "Please do your work here : " << __FILE__ <<" line " << __LINE__ <<std::endl;

  // The result is returned to the user
  return C;
}

/*
dmatrix_denseCM mulV2(const dmatrix_denseCM &A, const dmatrix_denseCM &B){
   // get the size data from A and B
  const size_t M = A.getNbLines();
  const size_t K = A.getNbColumns();
  const size_t N = B.getNbColumns();
  if (K != B.getNbLines()){
    // The code will stop if there is a missmatch on the size of A and B
      std::cout <<"error can't compute product of 2 matrices : dimension missmatch." << std::endl;
      throw;
  }
  // The size are correct, the matrix C that will store the result is created
  dmatrix_denseCM C(M,N,0.);
  // This is where you need to work : implement your matrix multiplication below, before return C

  //std::cout << "Please do your work here : " << __FILE__ <<" line " << __LINE__ <<std::endl;
  for(int j = 0; j < N; ++j)
  {
  	for(int k = 0; k < K; ++k)
  	{
  		for(int i = 0; i < M; ++i)
  		{
  			C(i,j) += A(i,k)*B(k,j);
  		}
  	}
  }
  // The result is returned to the user
  return C;
}
*/

dmatrix_denseCM mulV2(const dmatrix_denseCM &A, const dmatrix_denseCM &B)
{
	// get the size data from A and B
	const size_t M = A.getNbLines();
	const size_t K = A.getNbColumns();
	const size_t N = B.getNbColumns();
	if (K != B.getNbLines())
	{
		// The code will stop if there is a missmatch on the size of A and B
		std::cout <<"error can't compute product of 2 matrices : dimension missmatch." << std::endl;
		throw;
	}
	// The size are correct, the matrix C that will store the result is created
	dmatrix_denseCM C(M,N,0.);

	const double * a = A.data();
	const double * b = B.data();
	double * c = C.data();

	// This is where you need to work : implement your matrix multiplication below, before return C

	for(int k = 0; k < K; ++k)
	{
		for(int j = 0; j < N; ++j)
		{
			for(int i = 0; i < M; ++i)
			{
				*(c+i+j*M) += (*(a+i+k*M)) * (*(b+k+j*K));
			}
		}
	}

	//std::cout << "Please do your work here : " << __FILE__ <<" line " << __LINE__ <<std::endl;

	// The result is returned to the user
	return C;
}

dmatrix_denseCM muldgemm(const dmatrix_denseCM &A, const dmatrix_denseCM &B ){
   // get the size data from A and B
  const size_t M = A.getNbLines();
  const size_t K = A.getNbColumns();
  const size_t N = B.getNbColumns();
  if (K != B.getNbLines()){
    // The code will stop if there is a missmatch on the size of A and B
      std::cout <<"error can't compute product of 2 matrices : dimension missmatch." << std::endl;
      throw;
  }
  // The size are correct, the matrix C that will store the result is created
  dmatrix_denseCM C(M,N,0.);
  // This is where you need to work : implement your matrix multiplication below, before return C
  
  
  dgemm_( 'N', 'N', M, N, K, 1, A.data(), M, B.data(), K, 0, C.data(), M);
  
  //std::cout << "Please do your work here : " << __FILE__ <<" line " << __LINE__ <<std::endl;

  // The result is returned to the user
  return C;
}

double norm2(const dmatrix_denseCM &A ){
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

