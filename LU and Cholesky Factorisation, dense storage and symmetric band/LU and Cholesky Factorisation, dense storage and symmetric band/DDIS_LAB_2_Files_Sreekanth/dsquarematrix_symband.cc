#include "dsquarematrix_symband.h"
#include <cmath>

/**************************************!
   YOU HAVE TO IMPLEMENT THE FOLLOWING FUNCTION !!!
   !**************************************/  

double & dsquarematrix_symband::operator()(int i, int j) 
{
	if(i == j)
	{
		return *(a+((lb+1)*(i+1) -1));
	}
	else if (std::fabs(i-j) < lb)
	{
		int x = std::fabs(i-j) ;
		return  *(a+(lb+1)*std::max(i,j)+lb-x);
	}
	else
	{
		static double default_value = 0.0;
		return default_value;;
	}
}


dsquarematrix_symband::dsquarematrix_symband():a(nullptr ), m(0), lb(0){}

dsquarematrix_symband::dsquarematrix_symband(const dsquarematrix_symband &Ab):a(new double[Ab.m*(Ab.lb+1)] ), m( Ab.m), lb(Ab.lb){
    std::copy(Ab.a, Ab.a+m*(lb+1), a );
  }

dsquarematrix_symband::dsquarematrix_symband( const dmatrix_denseCM & A, int lb):a(nullptr),m(0), lb(0)
{
    set_size(A.getNbLines(), lb);
    std::fill(a,a+m*(lb+1),0.);
    for (int idiag = 0; idiag <=lb; ++idiag){
      for (int j = idiag; j< m; ++j){
	int i = j-idiag;
	(*this)(i,j) = A(i,j);
      }
    }
 }

void dsquarematrix_symband::set_size(int _m , int _lb){
	m = _m;
	lb = _lb;
	if (!a) delete []a;
	a = new double[(lb+1)*m];
}

 
  
double dsquarematrix_symband::operator()(int i, int j) const{
  if (abs(j-i) > lb) return 0;
  double res = (const_cast< dsquarematrix_symband &>( (*this) ))(i,j);
  return res;
}

dsquarematrix_symband::~dsquarematrix_symband(){
  delete []a;
}

const double * dsquarematrix_symband::data() const{
  return a;
}

double * dsquarematrix_symband::data() {
  return a;
}

int dsquarematrix_symband::getNbLines()  const {return m;}
int dsquarematrix_symband::getBandwidth() const {return lb;}
const double *dsquarematrix_symband::begin() const{return a;}
const double *dsquarematrix_symband::end() const{return a+m*(lb+1);}
double *dsquarematrix_symband::begin() {return a;}
double *dsquarematrix_symband::end() {return a+m*(lb+1);}
//output 
std::ostream& operator<<(std::ostream& out, const dsquarematrix_symband& A)
{
  const int M = A.getNbLines();
  const int lbu = A.getBandwidth();
  std::cout << M << ", " << M << " square symetric in compressed band storage, with lbu = "<< lbu <<std::endl;
  for (int ib = 0; ib < lbu+1; ++ib){
    for (int j = 0; j < M; ++j){
      out << A.data()[ib+j*(lbu+1)] << " ";
    }
    out << std::endl;
  }
  return out;
}


