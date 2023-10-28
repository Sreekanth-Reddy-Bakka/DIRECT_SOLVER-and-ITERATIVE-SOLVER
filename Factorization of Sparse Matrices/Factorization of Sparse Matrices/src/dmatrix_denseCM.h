#ifndef _dmatrix_denseCM_
#define _dmatrix_denseCM_
#include <iostream>
#include <cassert>
#include <algorithm>
#include <vector>

enum storage {CM, RM};

/// This class implement m*n matrix off double precision float in column major storage. Indexing start at 0.
class dmatrix_denseCM{
public:
  /// Construct an m*n matrix of double. All the terms are set to val.
  dmatrix_denseCM(size_t _m, size_t _n, double val);
  /// Construct an m*n matrix of double. No initialisation of the terms of the matrix
  dmatrix_denseCM(size_t _m=0 , size_t _n=0);
  /// This will construct a matrix of size mn and fil it with the content stored in the 1d array pointed to by a. if s = RM, the data in a are interpreted as being stored row by row if s =CM they are interpreted column by column 
  dmatrix_denseCM(size_t _m, size_t _n, const double *a, storage s=RM);
   /// Copy Constructor. typical usage : if A is a dmatrix_denseCM, dmatrix_denseCM B(A); will create B as a deep copy of A.
  dmatrix_denseCM(const dmatrix_denseCM &in);
  /// Assignment operator. typical usage : A = B;
  /*! A can be reallocate */
  dmatrix_denseCM& operator = (const dmatrix_denseCM &in);
  /// dmatrix_denseCM destructor. called automatically when the object goes out off scope
  ~dmatrix_denseCM();
  /// give read/write access to term i,j.
  /*!  usage : if A is a matrix (non const), A(2,3) = 2. set the term a_(2,3) to 2.; 
       double a =  A(2,3); set a to the value of a_(2,3).
       Throw an exeption if i>=m or j >=n, unless flag NDEBUG is set at compile time
  */
  double &operator()(size_t i, size_t j);
   /// give read access to term i,j.
  /*!  usage : if A is a matrix;
       double a =  A(2,3); set a to the value of a_(2,3).
       Throw an exeption if i>=m or j >=n, unless flag NDEBUG is set at compile time
  */
  const double & operator()(size_t i, size_t j) const;
  /// return the nuber of lines in a matrix. typical usage : m = A.getNbLines()  
  int getNbLines() const;
  /// return the nuber of lines in a matrix. typical usage : n = A.getNbColumns()   
  int getNbColumns() const;
  /// return a pointer to the beginning of the raw data of the matrix.
  /*! Note that this address become invalid it the matrix is reallocated for some reason 
    This is the "Const" version, the data can not be changed using the returned pointer.
  */
  const double * data() const;
  /// return a pointer to the beginning of the raw data of the matrix.
  /*! Note that this address become invalid it the matrix is reallocated for some reason 
     This is the no "Const" version: the data can  be changed using the returned pointer.
  */
  double * data();
  /// same as double *data() . here for compatibility reason to be more "stl" like
  double *begin();
  /// Return a pointer to one past the last term in the matrix
  double *end();
  const double *begin() const;
  const double *end()const;
  /// Add to the current matrix the matrix on the right hand side : A +=B
  /*! Throw an exeption if A and B are not the same size */
  dmatrix_denseCM &operator+=(const dmatrix_denseCM &r);
  /// Remove to the current matrix the matrix on the right hand side : A -=B
  /*! Throw an exeption if A and B are not the same size */
  dmatrix_denseCM &operator-=(const dmatrix_denseCM &r);
  /// Multiply all the term of the matrix by a : A *=a
  dmatrix_denseCM &operator*=( double a);
  /// Divide all the term of the matrix by a : A /= a
  dmatrix_denseCM &operator/=( double a);
  /// return sqrt (sum a_ij^2)
  double norm2() const;
  /// Reset the matrix to a  new size.
  void reallocate(int m, int n);
  /// set data using the value between begin and end
  /*! will throw if std::distance(begin,end)>= m*n 
      data between begin and end are considered as a collection of rows if s =RM
      or a collection of columns if s = CM
  */ 
  void setData(const double *begin, const double *end, storage s=RM );
 private:
  size_t m,n;
  double *a;

};

/// generate a random matrix of size m,n, with terms betwin min  and max
dmatrix_denseCM generateRandom(int m, int n, double min, double max);

/// matrix difference usage A = B-C
/*! throw if B and C are not the same size */
dmatrix_denseCM operator-(const dmatrix_denseCM &l, const dmatrix_denseCM &r );

/// matrix addition. usage A = B + C 
/*! throw if B and C are not the same size */
dmatrix_denseCM operator+(const dmatrix_denseCM &l, const dmatrix_denseCM &r );

/// matrix multiplication with a scalar  A = B *a;
dmatrix_denseCM operator*( const dmatrix_denseCM &l, double a );
/// matrix multiplication with a scalar  A = a *B;
dmatrix_denseCM operator*( double a, const dmatrix_denseCM &l );

/// ijk version of the product A*B
dmatrix_denseCM mulijk(const dmatrix_denseCM &A, const dmatrix_denseCM &B);

/// jki version of the product A*B
dmatrix_denseCM muljki(const dmatrix_denseCM &A, const dmatrix_denseCM &B);

/// product using dgemm
dmatrix_denseCM muldgemm(const dmatrix_denseCM &A, const dmatrix_denseCM &B);

/// product A*B (use muldgemm)
dmatrix_denseCM operator*(const dmatrix_denseCM &A, const dmatrix_denseCM &B);

/// return the transpose of A
dmatrix_denseCM transpose(const dmatrix_denseCM &A);

/// permut lines and columns of A using the permutation given by p
dmatrix_denseCM sympermut(const dmatrix_denseCM &A, const std::vector<int > &p);

/// permut lines of A using the permutation given by p
dmatrix_denseCM permutlines( const dmatrix_denseCM &A, const std::vector<int > &p);

/// permut columns of A using the permutation given by p
dmatrix_denseCM ipermutlines( const dmatrix_denseCM &A, const std::vector<int > &p);



/// return the norm2 of matrix A : a = norm2(A)
double norm2(const dmatrix_denseCM &A );

/// push the data in the matrix A on the screen or in a file.
/*! Typical usage :
    to print on screeen :
    std::cout << A; 
    to write in a file :
    std::ofstream f("myfile"); f << A; */ 
std::ostream& operator<<(std::ostream& os, const dmatrix_denseCM& A);


/// Read data from a file or from the keyboard. A is usually reallocated.
/*! Typical usage
   to read fron screen :
   std::cin >> A;
   to read from a file, assuming the file data.mat exist :
   std::ifstream f("data.mat"); f >> A; */
   std::istream& operator>>(std::istream& in, dmatrix_denseCM& A);


/// Read  data from the stream in, supposed in MatrixMarketFormat.
void readMatrixMarket(std::istream &in, dmatrix_denseCM &A);
/// Read data from stream, supposed in home format (Row by row)
void readMyFormat(std::istream &in, dmatrix_denseCM &A);


/// You must implement this function : it must return the upper bandwidth of A
int computeBandwidthUp(const dmatrix_denseCM &A );

/// You must implement this function : it must return the lower bandwidth of A
int  computeBandwidthDown(const dmatrix_denseCM &A );

/// compute the number of Nonzero in A
inline int  computeNNZ(const dmatrix_denseCM &A ){
  int nnz = 0;
  std::for_each(A.begin(), A.end() , [&nnz](double a){ if (a !=0.)  ++nnz;});
  return nnz;
}

/// return the number of non zero in A if a was stored in skyline format.
inline int  skylineNNZ(const dmatrix_denseCM &A ){
  int nnz = 0;
  int m = A.getNbLines();
  int n = A.getNbColumns();
  
  for (int l =0; l < m; ++l){
    int lnnz = m-1;
    while ((A(l,lnnz) ==0.)) --lnnz;
    nnz += lnnz;
  }
  return nnz;
}






#endif


