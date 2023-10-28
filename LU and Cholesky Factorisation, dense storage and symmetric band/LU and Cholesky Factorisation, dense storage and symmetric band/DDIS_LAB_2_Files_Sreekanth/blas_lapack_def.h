#ifndef _blas_lapack_def_
#define _blas_lapack_def_
extern "C"{
  // blas
  //L1
  void dswap_(const int &N, const double * DX, const int &INCX, const double * DY, const int &INCY );
  void dscal_(const int & N, const double & ALPHA, double * X, const int & INCX);
  //L2
  void dger_( const int &M, const int &N, const double &alpha,  const double *X, const int &INCX, double*Y, const int &INCY,  const double *A, const int &LDA);
  //l3
  void dgemm_( const char &TRANSA, const char &TRANSB, const int &M, const int  &N, const int & K,  const double & ALPHA,  const double *A, const int & LDA, const double *B,  const int &LDB, const double &BETA, double *C, const int & LDC);

  void dtrsm_(const char & SIDE, const char &UPLO, const char &TRANSA, const char &DIAG, const int &M, const int &N, const double &alpha, const double *A, const int &LDA, const double*B, const int &LDB);

  // lapack
  void dgetrf_( const int &M, const int &N,  double *A, const int &LDA, int *IPIV, int &INFO);
  void dgetrf2_( const int &M, const int &N,  double *A, const int &LDA, int *IPIV, int &INFO);
  void dgetf2_( const int &M, const int &N,  double *A, const int &LDA, int *IPIV, int &INFO);
  void dpbtrf_(const char &UPLO, const int &N, const int &KD, double *AB,  const int &LDAB, int &INFO);
  void dpbtrs_(const char &UPLO, const int &N, const int &KD, const int &NRHS, double *AB,  const int &LDAB, double *B, const int & LDB,  int &INFO);
  
}
#endif
