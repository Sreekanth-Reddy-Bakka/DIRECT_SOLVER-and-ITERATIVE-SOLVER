

#include "dmatrix_CCS.h"
#include "graph.h"
#include "mmio.h"

#include <cassert>

dmatrix_CCS::dmatrix_CCS( int _m, int _n, int _nnz, int *_lineindex, int *_columnptr, double *_a):
  m(_m), n(_n), nnz(_nnz), lineindex(nnz), columnptr(n+1), a(nnz){
  std::copy(_lineindex, _lineindex+nnz, lineindex.begin());
  std::copy(_columnptr,  _columnptr+n+1, columnptr.begin());
    if (_a==nullptr) std::fill(a.begin(), a.end(), 0.);
    else  std::copy(_a, _a+nnz, a.begin());
}

dmatrix_CCS::dmatrix_CCS( const graph & columng, double val  ){
    compressGraph( columng, columnptr, lineindex );
    m = n = columnptr.size()-1;
    nnz = columnptr[m];
    a.resize(nnz);
    std::fill(a.begin(), a.end(), val);
  }


std::istream& operator>>(std::istream& in, dmatrix_CCS& A){
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
	  printf(" Can't read the M, N and nz");
	  exit(1);
	}
	
	
	struct ija{
	  int i; int j; double a;
	};

	std::vector<ija> acoo;
	
	for (int ke =0; ke < nz; ++ke){
	  int i, j;
	  double v;
	  in >> i >> j >> v;
	  ija aij  {i-1,j-1,v };
	  acoo.push_back(aij);
	  if ( (aij.i != aij.j ) && (mm_is_symmetric(matcode))){
	    ija aijs {j-1, i-1, v};
	    acoo.push_back(aijs);
	  }
	}
	int nnz = acoo.size(); 
	auto lessthan =  [] (const ija &a, const ija &b){ 
		    if (a.j < b.j) return true;
		    if (a.j > b.j) return false;
		    if (a.i < b.i) return true;
		    return false;
	};
	auto equal = [&lessthan]( const ija &a, const ija& b){
	  return ( (!lessthan(a,b)) && (!lessthan(b,a)));
	};
	std::sort(acoo.begin(), acoo.end(), lessthan);
		 
	
	// Should check for multiple ij entry
        // Not now ...
	//int 
	acoo.erase(std::unique(acoo.begin(), acoo.end(), equal), acoo.end());
	if (nnz != acoo.size()){
	  std::cout << " WARNING : The input matrix have some term ij defined a multiple times. CHECK YOUR DATA   "<< std::endl;
	  throw;
	}
	std::vector<double >  a(nnz,0.);
	std::vector<int >  columnptr(N+1, 0);
	std::vector<int >  lineindex(nnz, 0);
	columnptr[0] = 0;
	int j =0;
        int k = 0;
	for(auto aij : acoo){
	  if (j != aij.j) {
	    while(j!= aij.j){ ++j; columnptr[j]=k; 
	    }
	  }
	  a[k] = aij.a;
	  lineindex[k] = aij.i;
	  ++k;
	}
	while(j<N){++j; columnptr[j] =k; }
        A =dmatrix_CCS(M,N, nnz,  lineindex.data(),columnptr.data(),  a.data());
        return in;
      }
   
  }
  else throw;
}

std::ostream& operator<<(std::ostream& out, const dmatrix_CCS& A){
  if(&out ==&std::cout){
    out << A.m << " " << A.n << " " << A.nnz << std::endl;
    out << "columnptr :" <<std::endl;
    for (int k = 0; k < A.n+1; ++k) out << A.columnptr[k] << " ";
    out << std::endl;
    out << "lineindex :" <<std::endl;
    for (int k = 0; k < A.nnz; ++k) out << A.lineindex[k] << " ";
    out << std::endl;
    out << "a : " << std::endl;
    for (int k = 0; k < A.nnz; ++k) out << A.a[k] << " ";
    out << std::endl;
  }
  else{
    //assuming matrix market format
    out << "%%MatrixMarket matrix coordinate real general" << std::endl;
    out << A.m << " "<< A.n << " " << A.nnz << std::endl;
    for(int j =0; j < A.n; ++j)
      for (int k =  A.columnptr[j];  k < A.columnptr[j+1]; ++k)
	out << A.lineindex[k]+1 << " " << j+1 << " "<< A.a[k] << std::endl; 
  }
  return out;
}
/*
  void getCol(int j, int &size, double *colb, int *lineid ){
  int scolj = colptr[j];
  int ecolj = colptr[j+1];
  size = endcolj - startcolj;
  colb = &a[scolj];
  lineid = &lineindex[scolj];
  }*/

/// will throw if trying to modify an a structurally zero term.
double &  dmatrix_CCS::operator()(int i, int j){
  assert( (i <m) && (j <n));
  int sj = columnptr[j];
  int ej = columnptr[j+1];
  auto it = std::find( &lineindex[sj], &lineindex[ej], i );
    if (it== &lineindex[ej]) throw;
    return a[ std::distance(&lineindex[0], it ) ];
}

double   dmatrix_CCS::operator()(int i, int j) const{
  int sj = columnptr[j];
  int ej = columnptr[j+1];
  auto it = std::find( &lineindex[sj], &lineindex[ej], i );
  if (it== &lineindex[ej]) return 0.;
  return a[ std::distance(&lineindex[0], it ) ];
}

int   dmatrix_CCS::getNbLines()const   {return m;}

int   dmatrix_CCS::getNbColumns()const {return n;}
  
/// destructor
/// Allocate: clear whatever was previously stored in the matrix, and then allocate storage for an mxn matrix with nnz non zero terms. owner is set to true : the user do not have to take care of the destruction of the newed array.
//void   dmatrix_CCS::allocate(const int& _m, const int & _n, const int &_nnz);
/// clear : matrix is set to empty. if owner, all array are destroyed

void   dmatrix_CCS::clear(){
  m = n = nnz =0;
  a.clear();
  columnptr.clear();
  lineindex.clear();
};


dmatrix_denseCM operator*( const dmatrix_CCS &A, const dmatrix_CCS &B){
  int M = A.getNbLines();
  int N = B.getNbColumns();
  int K = B.getNbLines();
  assert(K == A.getNbColumns());
  dmatrix_denseCM C( A.getNbLines(), B.getNbColumns(),0.);
  for (int k = 0; k < N; ++k){
    int sk = A.columnptr[k];
    int ek = A.columnptr[k+1];
    for (int kk = sk; kk < ek; ++kk){
      int i = A.lineindex[kk];
      double aik = A.a[kk];
      for (int j = 0; j < N; ++j){
	C(i,j) += aik*B(k,j);
      }
    }
  }
  return C;
}
  
dmatrix_denseCM trsm( const dmatrix_CCS &L, const dmatrix_denseCM &B){
  dmatrix_denseCM X = B;
  int N = L.getNbLines();
  int M = B.getNbColumns();
  for (int k = 0; k < M; ++k){ 
    for (int j = 0; j < N; ++j){
      X(j,k) = X(j,k)/L(j,j);
      for (int i = j+1; i< N; ++i ){
	X(i,k) -= L(i,j)*X(j,k);
      }
    }
  }
  return X;
}
/*
dmatrix_CSS generateRandomDominanteDiag_dmatrix_CSS(int m, int averagennzpercol){
  std::random_device rd;
  std::mt19937 gen(rd());  // this is the generator
  std::vector<int > columnindex(m+1);
  std::uniform_int_distribution<> dis(1, averagennzpercol*2+1);
  
  columnptr[0] = 0;
  int nnz = 0;
  for (int c =1; c < m+1; ++c){
    int n = dis(gen);
    nnz +=n;
    columnptr[c] = columnptr[c-1] + n;
  }
  
  std::uniform_double_distribution<> dis_a(0.,1);
  std::uniform_int_distribution<>    dis_l(0, m);
  for (int c =0; c < m; ++c){
    for( int ii = 0; ii < columnptr[c+1] - columnptr[c]; ++ii){
      a[columnindex[c + ii] ] = dis_a(gen);
      lineindex[columnindex[c + ii] ] = dis_a(gen);
    }
    
  }
  
}

dmatrix_CSS generateRandomLowerTriangularUnitDiag_dmatrix_CSS(int m, int averagennz){
  std::random_device rd;
  std::mt19937 gen(rd());  // this is the generator

  std::vector<int > columnindex(m+1);
  columnindex[0] = 0;
  std::uniform_int_distribution<> dis(1, m-1);
  

  

}
*/


dmatrix_denseCM convertToDenseCM(const dmatrix_CCS & in){
  int M = in.getNbLines();
  int N = in.getNbColumns();
  dmatrix_denseCM A(M,N,0.);
  for (int j =0; j < N; ++j){
    int sj = in.columnptr[j];
    int ej = in.columnptr[j+1];
    for (int kk = sj; kk < ej; ++kk){
      int i = in.lineindex[kk];
      A(i,j) = in.a[kk];
    }
  }
  return A;
}

dmatrix_CCS convertToCCS(const dmatrix_denseCM &in, double eps){
  int M = in.getNbLines();
  int N = in.getNbColumns();
  assert(M == N);
  auto gA = buildColumnGraph(in, eps);
  assert( gA.nbNodes() == N);
  int NNZ = gA.nbEdges();
  dmatrix_CCS A(M,N,NNZ);
  compressGraph(gA, A.columnptr, A.lineindex);
  
  for (int j =0; j < N; ++j){
    int sj = A.columnptr[j];
    int ej = A.columnptr[j+1];
    for (int kk = sj; kk < ej; ++kk){
      int i = A.lineindex[kk];
      A.a[kk] = in(i,j);
    }
  }
  
  return A;
}

dmatrix_CCS extractUpperTriangle(const dmatrix_CCS &A, int rmdiag){
  int N = A.n;
  int M = A.m;
  assert (N== M);
  std::vector<int> nzcol(N); // numper of nz in the colomn j, below diag discarded.
  for (int j= 0; j < N; ++j){
    for(int k = A.columnptr[j]; k < A.columnptr[j+1]; ++k ){
      int i = A.lineindex[k];
      if ( i<= j-rmdiag) nzcol[j]++; 
    }
  }
  /* std::cout << "NZCOL " << std::endl;
  for( int nz : nzcol){
    std::cout << nz << " ";
  }
  std::cout << std::endl;
  */
  std::vector<int > columnptr(N+1,0);
  columnptr[0] = 0;
  
  // c++17 scan ..
  for(int k=0; k< N; ++k)  columnptr[k+1] = columnptr[k]+nzcol[k];
  int NNZ = columnptr[N];
  std::vector<int > lineindex(NNZ,0);
  std::vector<double > a(NNZ,0.);
  int kk = 0;
  for (int j= 0; j < N; ++j){
    for(int k = A.columnptr[j]; k < A.columnptr[j+1]; ++k ){
      int i = A.lineindex[k];
      if ( i<= j-rmdiag){
	a[kk] = A.a[k];
	lineindex[kk] = i;
	++kk;
      }  
    }
  }
  return dmatrix_CCS(M,N, NNZ, lineindex.data(), columnptr.data(), a.data() );
}

int lbdown ( const dmatrix_CCS & A){
  int lb = 0;
  int N = A.n;
  for (int j = 0; j < N ;++j){
    for(int k = A.columnptr[j]; k < A.columnptr[j+1]; ++k ){
      int i = A.lineindex[k];
      if (i> j) lb = std::max(lb,(i-j));
    }
  }
  return lb;
}


dmatrix_CCS sympermut(const dmatrix_CCS &A, const std::vector<int > &p){
  int m = A.m;
  assert( (m == A.n) && (m == p.size()));
  std::vector<int > pt(m);
  for (int i =0; i < m; ++i) pt[ p[i] ] = i;
  int nnz = A.nnz;
  dmatrix_CCS pAp(m,m, nnz);
  pAp.columnptr[0] = 0;
  for (int j = 0; j< m ;++j){
    int aj = p[j];
    int baj = A.columnptr[aj];
    int eaj = A.columnptr[aj+1];
    int bpaj = pAp.columnptr[j];
    int epaj = bpaj+eaj-baj;
    pAp.columnptr[j+1] = epaj;
    int kpa = bpaj;
    for (int ka = baj; ka < eaj; ++ka,++kpa){
      int ia =  A.lineindex[ka];
      int ipa = pt[ia];
      pAp.lineindex[kpa] = ipa;
    }
    std::sort(&pAp.lineindex[bpaj], &pAp.lineindex[epaj]);
  }
  
  for(int j = 0; j < m; ++j){
    for(int k = pAp.columnptr[j]; k < pAp.columnptr[j+1]; ++k){
      int i = pAp.lineindex[k];
      int Ai   = p[i];
      int Aj   = p[j]; 
      pAp(i,j) = A(Ai,Aj);
    }
  }
  return pAp;
}

dmatrix_CCS transpose(const dmatrix_CCS & A){
  int m = A.m;
  int n = A.n;
  int nnz = A.nnz;
  std::vector<int> rowcount(m, 0); //store thr number of non_zero in row i of A;
  for(int k =0; k <nnz; ++k) {rowcount[ A.lineindex[k]]++;}
  std::vector<int> rowpointer(m+1,0);
  for(int i=0; i < m; ++i) rowpointer[i+1] = rowpointer[i]+rowcount[i];
  std::vector<int > nextinrow(rowpointer);
  std::vector<int > lineindex(nnz);
  std::vector<double > at(nnz);
  for (int j = 0 ;  j < n ; j++){
    for (int p = A.columnptr [j] ; p < A.columnptr [j+1] ; p++){
      int Ai = A.lineindex[p];
      int q = nextinrow[Ai]++;
      lineindex[q] = j;
      at[q] = A.a[p];
    }
  }
  return dmatrix_CCS(m,n, nnz, lineindex.data(), rowpointer.data(), at.data()  );
}



dmatrix_denseCM solveL ( const dmatrix_CCS &L, const dmatrix_denseCM &B ){
  int n = L.n;
  int LDB =B.getNbLines(); 
  assert((B.getNbLines() == n));
  int nrhs = B.getNbColumns();
  dmatrix_denseCM X(B);
  double *x = X.data();
  const double *Lx = L.a.data();
  const int    *Lp = L.columnptr.data();
  const int    *Li = L.lineindex.data();
  
  for (int j = 0 ; j < n ; j++){
    const double tmp = Lx [Lp [j]] ;
    for( int  k = 0; k <  nrhs; ++k ) x [j+k*LDB] /=tmp ;
    for (int p = Lp [j]+1 ; p < Lp [j+1] ; p++){
      const int Lip = Li[p];
      const double Lxp = Lx[p];
      for(int k = 0; k <nrhs; ++k ) x [Lip + k*LDB] -= Lxp * x [j + k*LDB] ;
    }
  }
  return X;
}





dmatrix_denseCM solveLT ( const dmatrix_CCS &L, const dmatrix_denseCM &B ){
  int n = L.n;
  int LDB =B.getNbLines(); 
  assert((B.getNbLines() == n));
  int nrhs = B.getNbColumns();
  dmatrix_denseCM X(B);
  double *x = X.data();
  const double *Lx = L.a.data();
  const int    *Lp = L.columnptr.data();
  const int    *Li = L.lineindex.data();
  
  for (int j = n-1 ; j >= 0 ; --j){
    int p;
    for ( p = Lp [j]+1 ; p < Lp [j+1] ; ++p){
      const int Lip = Li[p];
      const double Lxp = Lx[p];
      for(int k = 0; k <nrhs; ++k ) x [j + k*LDB] -= Lxp * x [Lip + k*LDB] ;
    }
    const double tmp = Lx [Lp [j]] ;
    for( int  k = 0; k <  nrhs; ++k ) x [j+k*LDB] /=tmp ;
  }
  return X;
}

dmatrix_denseCM solveU ( const dmatrix_CCS &U, const dmatrix_denseCM &B ){
  int n = U.n;
  int LDB = B.getNbLines();
  int nrhs = B.getNbColumns();
  assert((B.getNbLines() == n));
  
  dmatrix_denseCM X(B);
  double *x = X.data();
  const double *Ux = U.a.data();
  const int    *Up = U.columnptr.data();
  const int    *Ui = U.lineindex.data();
  for (int j = n-1 ; j >= 0 ; j--){
    double tmp =  Ux [Up [j+1]-1] ;
    for (int k = 0; k < nrhs; ++k) x [j+ LDB*k] /= tmp;
    for (int p = Up [j] ; p < Up [j+1]-1 ; p++){
      const int Uip = Ui[p];
      const double Uxp = Ux[p]; 
      for (int k = 0; k < nrhs; ++k) x[Uip + k *LDB] -= Uxp*x[j + k*LDB]; 
      //      x [Ui [p]] -= Ux [p] * x [j] ;
    }
  }
  return X;
}


dmatrix_denseCM solveUT ( const dmatrix_CCS &U, const dmatrix_denseCM &B ){
  int n = U.n;
  int LDB = B.getNbLines();
  assert((B.getNbLines() == n));
  int nrhs = B.getNbColumns();
  dmatrix_denseCM X(B);
  double *x = X.data();
  const double *Ux = U.a.data();
  const int    *Up = U.columnptr.data();
  const int    *Ui = U.lineindex.data();
  
  for ( int j = 0 ; j < n ; j++) {
    int p;
    for (  p = Up [j] ; p < Up [j+1]-1 ; p++){
      const double Uxp = Ux[p];
      const int Uip = Ui[p];
      for (int k = 0; k < nrhs; ++k) x[j+k*LDB] -= Ux [p] * x [Ui [p] + k*LDB] ;
    }
    const double Uxp = Ux[p];
    for (int k = 0; k < nrhs; ++k) x[j+k*LDB] /= Uxp;
  }
  
  return X;
}




