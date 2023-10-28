#include <iostream>
#include "dmatrix_CCS.h"
#include "slu_ddefs.h" //include file of the superLU library.

// This program gives you an exemple of using superLU

// you can find the superLU library and the (mostly ) up to date documentation here:
//http://crd-legacy.lbl.gov/~xiaoye/SuperLU/



int main(){
  // small test with a matrix stored in CCS.
  // solving Ax =b with superLU
  /* A = 1. 2. 0                
         0  3. 0                
         4. 0  5.        
     m = 3, n =3, nnz = 5
     b = {1, 1, 1}
     solution is :
     x = {1/3, 1/3, -2/3}
     
     CCS : compressed column storage.
     Accs : m =3, n =3, nnz = 5
     lineindex    = 0 2 0 1 2
     columnptr    = 0 2 4 5
     values       = 1. 4. 2. 3. 5.

  */
  
  int m = 3;
  int n = 3;
  int nnz = 5;
 

  int linei[]   = {0, 2, 0, 1, 2 };
  int colptr[]  = {0, 2, 4, 5};
  double a[]   =  {1.,4.,2.,3.,5.};
  dmatrix_CCS Accs(m,n,nnz);
    
  std::copy(a, a+nnz, Accs.a.begin() );
  std::copy(linei, linei+nnz, Accs.lineindex.begin() );
  std::copy(colptr, colptr+n+1, Accs.columnptr.begin() );
  
  // create a right hand side. here it's full of 1.
  int  nrhs   = 1;
  double * rhs = new double [m*nrhs];
  for (int i = 0; i < m ;++i) rhs[i] =1.;

  
  // To this point we have stored the matrix A in our own CCS format.
  
  // SuperLU use a design called "handle" to access to  matrix of diverse format. The handle to a matrix
  //  type is  called SuperMatrix. This type is defined by the superLU library.
  
  SuperMatrix A; // A a handle to the sparse matrix A
  SuperMatrix B; // B a handle to a dense matrix for the rights hand side, and the solution uppon completion
  
  SuperMatrix L; // L :  a handle to the Lower triangular part of the factorization. the solver will allocate it
  SuperMatrix U; // U :  a handle to the Upper triangular part of the factorization. the solver will allocate it
  
  // Fill up the handle A in super LU format. Note that all the relevent pointer are passed to the function.
  // The "type" of A is defined by the constants passed at the  end of the function : SLU_NC means compressed column, SLU_D means the matrix contains double, and SLU_GE means the matrix is General (in the sense that it is not in symmetric storage for instance)
  dCreate_CompCol_Matrix(&A, m, n, nnz, Accs.a.data(), Accs.lineindex.data(), Accs.columnptr.data(), SLU_NC, SLU_D, SLU_GE);


   // Fill up the handle B in super LU format. Note that all the relevent pointer are passed to the function.
  // The "type" of B is defined by the constants passed at the  end of the function : SLU_DN means dense matrix column major format, SLU_D means the matrix contains double, and SLU_GE means the matrix is General (in the sens that it is not symmetric for instance)
  dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
  
  // This variable store the different options for the superLU solver. Some options are only relevant for certains function of the superLU library.
  superlu_options_t options;  

  // Set the default options:
  set_default_options(&options);
  /*  The default option are :
		       	 options.Fact = DOFACT;
	       		 options.Equil = YES;
       			 options.ColPerm = COLAMD;
       			 options.DiagPivotThresh = 1.0;
       			 options.Trans = NOTRANS;
	       		 options.IterRefine = NOREFINE;
	       		 options.SymmetricMode = NO;
	       		 options.PivotGrowth = NO;
	       		 options.ConditionNumber = NO;
	       		 options.PrintStat = YES;
	        */
  /* Some explaination of some option below .. Up to you to modify them and see the effect
     options.Fact = DOFACT/ FACTORED/ SamePattern / SamePattern_SameRowPerm
           The Fact option refere to the work to be done by some of the "driver" function.
            DOFACT : DO the  factorization
            FACTORED : the factorizationis already done, use the factor as is for whatever need to be done (like solve ...)
	    SamePattern : The Matrix has the same non zero pattern as the previous one (typical of a non linear finite element or finite volume problem : we have to solve for different matrix that have the same patterns for different time step or inside iterations of a newton procedure.) The idea is to reuse the same symbolic factorization (this is only avalaible in expert routine)
	    SamePattern_SameRowPerm :Row Perm  is used to select best pivot (pivoting). If we know that the new matrix we want to solve for is similar to the previous one, we can use the same pivoting strategie.


    options.Equil = NO / YES; : Equil refere to the equilibration of the system. Lines of the matrix are multiplied by some value before the factorisation in order to improve on the quality of the solution, and limit round of error.
     options.ColPerm = NATURAL/ MMD_ATA / MMD_AT_PLUS_A /COLAMD ColPerm refer to the renumbering strategy used to reduce the Fill in. There is also the option to give your own renumbering, not listed here
     options.DiagPivotThresh = 0.0 ... 1.0 This parameter control the pivoting. if set to 0.0, almost no pivoting will be applied, 1.0, lot of pivoting.

    options.PivotGrowth = YES / NO compute the growth of pivot : usefull to asses quality of the solution;
    options.ConditionNumber = YES/ NO compute the condition number (an approximation of it using the LU factorization ); 
    options.IterRefine = NOREFINE / SLU_SINGLE / SLU_DOUBLE /SLU_EXTRA in expert mode, this option refere to iterative refinement of the solution, with different level of refinment. Iterative refinment of the solution iteratively improve the solution using an iterative solver viewing the LU factor as a way to approximate the inverse. It also turn on in expert mode, the computation of forward and backward error estimation of the solution.

  */

  
  SuperLUStat_t stat; // this variable is used to store some statistics.
  StatInit(&stat);
  int info; // this variable contains informations on the succes or failure of the superLU fonction
  int *perm_c = new int [n]; // This variable stores the permutations applied on the columns of the matrix
  int *perm_r =  new int[m]; // This Variable stores the permutations applied on the row of the matrix
  

  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);


  std::cout << "Solution :" << std::endl;
  for (int i=0; i < m; ++i)  std::cout << " "<<  rhs[i] ;
  std::cout << std::endl;
  
  //retreive some statistics on the factors
  NCformat       *Ustore;
  SCformat       *Lstore;
  Lstore = (SCformat *) L.Store;
  Ustore = (NCformat *) U.Store;
  printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
  printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
  printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
  printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);

  if ( options.PrintStat ) StatPrint(&stat);
  
  // Solving again using the same factor.
  dgstrs (NOTRANS,  &L, &U, perm_c, perm_r, &B, &stat, &info); 

  
  //cleaningup dynamically allocated variable.
  StatFree(&stat);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  Accs.clear();
  
  delete []perm_r;
  delete []perm_c;
  delete []rhs;

  
}
