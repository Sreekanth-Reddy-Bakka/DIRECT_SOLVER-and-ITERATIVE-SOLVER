#ifndef _dmatrix_CCS_
#define _dmatrix_CCS_

#include "dmatrix_denseCM.h"

class graph;

class dmatrix_CCS{
 public:
  /// m :  number of line, n : number of columns
  int m,n;
  /// nnz number of non zero terms
  int nnz;
 
  /// lineindex[k] contain the line number of the kth non zero term in a.
  std::vector<int > lineindex;
  /// columnptr[j] contain the index k of the beginning of column j in a and lineindex. columns
  std::vector<int >  columnptr;
  /// pointer to the beginning of the array containing the non zero term of the matrix.
  std::vector<double > a ; 
  /// owner is a bool that state if an object of the class own its memory (in this case, all the array will be deleted at destruction of an object of this class)
   inline dmatrix_CCS( ):m(0), n(0),nnz(0), lineindex(nnz,0), columnptr(n+1,0), a(nnz,0.){}
  inline dmatrix_CCS( int _m, int _n,  int _nnz):m(_m), n(_n),nnz(_nnz), lineindex(nnz,0), columnptr(n+1,0), a(nnz,0.){}
 
  dmatrix_CCS( const graph & columng, double val =0. );
  
  dmatrix_CCS(  int _m, int _n, int _nnz, int *_lineindex, int *_columnptr, double *_a= nullptr );
  
  
  /// will throw if trying to modify an a structurally zero term.
  double &operator()(int i, int j);
  
  /// return 0. if aij not stored.
  double operator()(int i, int j) const;

  /// return number of lines
  int getNbLines()const;
  
  /// return number of Columns
  int getNbColumns()const;
  
  
  
  /// Allocate: clear whatever was previously stored in the matrix, and then allocate storage for an mxn matrix with nnz non zero terms. owner is set to true : the user do not have to take care of the destruction of the newed array.
  void allocate(const int& _m, const int & _n, const int &_nnz);
  /// clear : matrix is set to empty. if owner, all array are destroyed
  void clear();
};

/// return the transpose of A in CCS format.
dmatrix_CCS transpose(const dmatrix_CCS & A);
/// return X solution of LX = B; where L is lower triangular
dmatrix_denseCM solveL ( const dmatrix_CCS &L, const dmatrix_denseCM &B );
/// return X solution of L'X = B; where L is lower triangular
dmatrix_denseCM solveLT ( const dmatrix_CCS &L, const dmatrix_denseCM &B );
/// return X solution of UX = B; where U upper triangular
dmatrix_denseCM solveU ( const dmatrix_CCS &U, const dmatrix_denseCM &B );
/// return X solution of U'X = B; Where U is upper triangular
dmatrix_denseCM solveUT ( const dmatrix_CCS &U, const dmatrix_denseCM &B );


/// Read the matrix from a file stream (.mtx format)
 std::istream& operator>>(std::istream& in, dmatrix_CCS& A);
 
/// write A on the screen (CCS style)
/*! To write on screen in a more "Readable" way
    you can convert first A in dense_CM format : 
     std::cout << convertToDenseCM(A) << std::endl;
 */
 std::ostream& operator<<(std::ostream& out, const dmatrix_CCS& A);


/// the following function would be use like: C = A*B, where C and B are dense and A is sparse. It is not implemented yet. That somethig for lab 4
//dmatrix_denseCM operator*( const dmatrix_CCS &A, const dmatrix_dense CM &B);

/// return a dense_CM matrix from the inputed CCS matrix
dmatrix_denseCM convertToDenseCM(const dmatrix_CCS & A);

/// return CCS matrix from the inputed denseCM matrix. Aij is considered null if | Aij | < eps 
dmatrix_CCS convertToCCS(const dmatrix_denseCM &A, double eps =0.);

/// permut matrix A by exchanging line and column according to the same permutatio vector p.
dmatrix_CCS sympermut(const dmatrix_CCS &A, const std::vector<int > &p);

/// extract the upperTriangle of of A. if rmdiag != 0 : rmdiag =1 remove 1 diag, rmdiag =m remvoe m diag rmdiag = -m add m subdiagonal
dmatrix_CCS extractUpperTriangle(const dmatrix_CCS &A, int rmdiag = 0);


/// return the lower band width
int lbdown ( const dmatrix_CCS & A);
#endif
