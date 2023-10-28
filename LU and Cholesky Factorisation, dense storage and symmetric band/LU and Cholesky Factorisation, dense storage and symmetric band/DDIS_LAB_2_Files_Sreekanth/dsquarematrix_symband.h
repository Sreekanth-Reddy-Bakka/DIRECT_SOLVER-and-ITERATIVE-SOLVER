#ifndef _dsquarematrix_symband_
#define _dsquarematrix_symband_
#include "dmatrix_denseCM.h"


/// This class is maint to store a symmetric matrix in band format
/*! It is quite minimimalistic since only the needed function for the exercise are defined 
About the assumed storage : compressed band symmetric Column major ...
For example, given the following symmetric 5x5 matrix :
50 2 3 0 0 
2 51 6 7 0 
3 6 52 9 8 
0 7 9 53 10 
0 0 8 10 54
m = 5 is the number oflines/columns.

The band width lb is 2, and it is stored Like a rectangular matrix of lb+1 = 3 lines and 5 columns, in a 1d array, following the column major convention:
x   x  3  7  8 
x   2  6  9 10 
50 51 52 53 54 
where the x represent useless value that are never addressed.
In CM storage, the 1d array looks like this :
x x 50 x 2 51 3 6 52 7 9 53 8 10 54, where the size of the 1d array is (lb+1)*m
*/
class dsquarematrix_symband{
 public :
  /// default constructor : create an empty matrix
  dsquarematrix_symband();
  /// copy constructor : make a new squarematrix_symband by copying Ab
  dsquarematrix_symband(const dsquarematrix_symband &Ab);
  /// Create a symband matrix, by taking it's term for the MAtrix A in denseCM storage. the number of line of the constructed matrix is the same as the one in A, lb is the bandwidth (number off non-zero superdiagonal considered.)
  /*   Only the term in the band are taken from A
     Note that it uses the operator()(i,j) of the band matrix. So it won't work properlly until you implemented this function 
 */ 
  dsquarematrix_symband( const dmatrix_denseCM & A, int lb);
  /// This function reallocate the matrix. _m is the number of line, and lb the bandwidth.
  void set_size(int _m , int _lb);  
  /// This return the i, j term of the matrix if it exist, else return 0. it can be called on a dsquarematrix_symband A by  calling A(i,j)
  /**************************************!
   YOU HAVE TO IMPLEMENT THIS FUNCTION IN THE 
    dsquarematrix_symband.cc  file.
   !**************************************/
  inline double &operator()(int i, int j);
  /// Const version of the above version. Implemented already in terms of the above function.  
  inline double operator()(int i, int j) const;
  /// destructor
  ~dsquarematrix_symband();
  /// return a pointer to the begining of the internal array
  const double * data() const;
  /// return a pointer to the begining of the internal array
  double * data();  
  /// return the number of lines (and columns)
  int getNbLines()  const;
  /// return the Bandwidth
  int getBandwidth() const;
  /// return a pointer to the begining of the internal array
  const double *begin() const;
  /// return a pointer one past the last term in the internal array
  const double *end() const;
  /// return a pointer to the begining of the internal array
  double *begin();
  /// return a pointer one past the last term in the internal array
  double *end();
 private :
  double *a;
  int m;
  int lb;
};
// output the matrix in band format
std::ostream& operator<<(std::ostream& os, const dsquarematrix_symband& A);
#endif
