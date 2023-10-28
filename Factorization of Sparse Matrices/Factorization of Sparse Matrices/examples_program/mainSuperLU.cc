
#include <iostream>
#include "slu_ddefs.h"
#include <fstream>
#include <chrono>
#include <mesh.h>
#include <random>

#include "dLU_denseCM.h"
#include "dLLT_denseCM.h"
#include "dmatrix_CCS.h"
#include "dLLT_CCS.h"
#include <queue>
#include "mat2gnuplot.h"
#include "dmatrix_denseCM.h"


dmatrix_CCS generateFromCommandLineArguments(int argc, char *argv[])
{
	std::string inputfile = "data/small_sym_sparse.mtx";
	bool printmat  = false;

	if (argc < 2){
		std::cout << "Usage : " << std::endl;
		std::cout << argv[0] << " input_filename" <<  " 0 " << std::endl;
		std::cout << "or :" << std::endl;
		std::cout << argv[0] << " input_filename" <<  " 1 " << std::endl;
		std::cout << std::endl<<  " Where input_filename is the input file, which can be a matrix format file (.mtx) or a gmsh mesh file (.msh ) " << std::endl;
		std::cout << " - If a mesh file is given, the analysis will be done on the classical finite element mass matrix build from the mesh " << std::endl;
		std::cout << " - If a mtx file is given, the analysis will be done directly on the matrix discribed by this file" << std::endl;
		std:: cout << " The second argument (0 or 1), is used to tell the code if you want to print the matrix shape on the screen (1), or not" << std::endl;
		std::cout << "Trying to run with default value : " << inputfile << " " << printmat << std::endl;
	}
	else{
		inputfile = argv[1];
		if (argc > 2)
			printmat = atoi(argv[2]);
	}
	bool from_mesh = (inputfile.npos != inputfile.find(".msh" ) );
	bool from_mtx = (inputfile.npos != inputfile.find(".mtx" ) );
		
	dmatrix_CCS A;
	if (from_mtx)
	{
		std::ifstream matfile(inputfile);
		if (!matfile.is_open()){ std::cout << "Can't open matrix file " << inputfile << std::endl; throw;}
		std::cout << "Reading Matrix from file " << inputfile << std::endl; 
		matfile >> A;
	}
	if (from_mesh)
	{
		std::ifstream meshin(inputfile);
		if (!meshin.is_open()) {std::cout << "Can't open mesh file " << inputfile << std::endl;   throw;}
		std::cout << "Reading Mesh from file " << inputfile << std::endl; 
		mesh m;
		meshin >> m;
		std::cout << m.nodes.size() << " "<< m.triangles.size() << std::endl;
		fem fem1(m);
		A = fem1.generateMassMatrix_CCS();
	}

	if (printmat) std::cout << "A : \n" <<  convertToDenseCM(A) << std::endl;
	std::cout << "matrix size : m = " << A.m<< " n = " << A.n << " nnz =" << A.nnz << std::endl;
	return A;
}


int main(int argc, char *argv[])
{
	/*
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.1, 1.);
	*/
	
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::time_point<std::chrono::system_clock> start2, end2;
	
  	dmatrix_CCS A = generateFromCommandLineArguments(argc, argv);
  	
	int m = A.getNbLines();
  	int n = A.getNbColumns();
  	int nnz = A.nnz;
  	
	//dmatrix_CCS A(A1.m, A1.n, A1.nnz);
	int nrhs = 3; 
	dmatrix_denseCM B = generateRandom(m,nrhs, 0., 1.);

/* ######################################################################################################
################                 CUTHILL MCKEE IN SUPERLU           #####################################
###################################################################################################### */
/*
	std::vector< int > neworder = cuthillmckee( buildSymGraph(A) );
	{ 
		std::cout << "Trying to solve the permuted system using cuthillmckee ordering " << std::endl;
		/// permuting the system  :  pA = PAPT, pB = PB
		dmatrix_CCS     pAp = sympermut(A, neworder);
		dmatrix_denseCM pB =  permutlines(B, neworder);
		//start = std::chrono::system_clock::now();

		dLLT_CCS LLT_pAp(pAp);
		
		//end = std::chrono::system_clock::now();
		



		//dmatrix_denseCM AT = convertToDenseCM(LLT_pAp.getLT());
		//dmatrix_denseCM A =  convertToDenseCM(LLT_pAp.getL());
		//dmatrix_denseCM B = A + AT;
		//spy_ge(B, "CUTHUILL_CUBE2.plt" );
		//spy_ge(B, "CUTHUILL_CUBE.plt" );
		//spy_ge(B, "CUTHUILL_sq.plt" );
		//spy_ge(B, "CUTHUILL_sq2.plt" );
		//spy_ge(B, "CUTHUILL_sq3.plt" );
		//spy_ge(B, "CUTHUILL_square4.plt" );
		//spy_ge(B, "CUTHUILL_b14.plt" );
		//spy_ge(B, "CUTHUILL_b15.plt" );

		//dmatrix_denseCM X  = ipermutlines(pX, neworder);
		SuperMatrix AA, BB, L, U;
		
		dCreate_CompCol_Matrix(&AA, pAp.m, pAp.n, pAp.nnz, pAp.a.data(), pAp.lineindex.data(), pAp.columnptr.data(), SLU_NC, SLU_D, SLU_GE);
			
		dCreate_Dense_Matrix(&BB, pAp.m, nrhs, pB.data(), m, SLU_DN, SLU_D, SLU_GE);
		
		superlu_options_t options;  

		// Set the default options:
		set_default_options(&options);
		

		SuperLUStat_t stat;             // this variable is used to store some statistics.
		StatInit(&stat);
		int info;                       // this variable contains informations on the succes or failure of the superLU fonction
		int *perm_c = new int [n];      // This variable stores the permutations applied on the columns of the matrix
		int *perm_r =  new int[m];      // This Variable stores the permutations applied on the row of the matrix


		 // Timing for factorization
		start2 = std::chrono::system_clock::now();
		dgssv(&options, &AA, perm_c, perm_r, &L, &U, &BB, &stat, &info);
		end2 = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds2 = end2-start2;
		std::cout << "Time taken to factorize using superLu with cuthill_reordering: " <<  (elapsed_seconds2).count() << std::endl;

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
		// Check for errors in factorization
		if (info == 0) {
		// Factorization successful
		// Timing for solving
			start = std::chrono::system_clock::now();
			dgstrs(NOTRANS, &L, &U, perm_c, perm_r, &BB, &stat, &info);
			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			//std::cout << "Time taken to solve Ax=B using SuperLU with cuthill re_ordering " << (elapsed_seconds).count() << std::endl;
		} else {
			std::cerr << "SuperLU factorization failed with error code: " << info << std::endl;
		}

		//cleaningup dynamically allocated variable.
		StatFree(&stat);
		Destroy_SuperMatrix_Store(&AA);
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
		Destroy_SuperMatrix_Store(&BB);
		A.clear();

		delete []perm_r;
		delete []perm_c;
		
	}
*/
/* ######################################################################################################
################                 REVERSE CUTHILL MCKEE IN SUPERLU           #####################################
###################################################################################################### */
/*
	std::vector< int > neworder = cuthillmckee( buildSymGraph(A) );
	{ 
		std::cout << "Trying to solve the permuted system using cuthillmckee ordering " << std::endl;
		/// permuting the system  :  pA = PAPT, pB = PB
		dmatrix_CCS     pAp = sympermut(A, neworder);
		dmatrix_denseCM pB =  permutlines(B, neworder);
	}
	{ 
		//std::cout << "Trying to solve the permuted system using reverse cuthillmckee ordering " << std::endl;
		std::vector< int > rorder(neworder.size());
		std::copy(neworder.rbegin(), neworder.rend(), rorder.begin());
		/// permuting the system  :  pA = PAPT, pB = PB
		dmatrix_CCS     pAp = sympermut(A, rorder);
		dmatrix_denseCM pB =  permutlines(B, rorder);
		start = std::chrono::system_clock::now();
		dLLT_CCS LLT_pAp(pAp);
		end = std::chrono::system_clock::now();         

		//dmatrix_denseCM AT = convertToDenseCM(LLT_pAp.getLT());
		//dmatrix_denseCM A =  convertToDenseCM(LLT_pAp.getL());
		//dmatrix_denseCM B = A + AT;
		//spy_ge(B, "ReverseCUTHUILL_CUBE2.plt" );
		//spy_ge(B, "ReverseCUTHUILL_CUBE.plt" );
		//spy_ge(B, "ReverseCUTHUILL_square.plt" );
		//spy_ge(B, "ReverseCUTHUILL_square2.plt" );
		//spy_ge(B, "ReverseCUTHUILL_square3.plt" );
		//spy_ge(B, "ReverseCUTHUILL_square4.plt" );

		//spy_ge(B, "ReverseCUTHUILL_b14.plt" );
		//spy_ge(B, "ReverseCUTHUILL_b15.plt" );
		SuperMatrix AA, BB, L, U;
		
		dCreate_CompCol_Matrix(&AA, pAp.m, pAp.n, pAp.nnz, pAp.a.data(), pAp.lineindex.data(), pAp.columnptr.data(), SLU_NC, SLU_D, SLU_GE);
			
		dCreate_Dense_Matrix(&BB, pAp.m, nrhs, pB.data(), m, SLU_DN, SLU_D, SLU_GE);
		
		superlu_options_t options;  

		// Set the default options:
		set_default_options(&options);
		

		SuperLUStat_t stat;             // this variable is used to store some statistics.
		StatInit(&stat);
		int info;                       // this variable contains informations on the succes or failure of the superLU fonction
		int *perm_c = new int [n];      // This variable stores the permutations applied on the columns of the matrix
		int *perm_r =  new int[m];      // This Variable stores the permutations applied on the row of the matrix


		 // Timing for factorization
		start2 = std::chrono::system_clock::now();
		dgssv(&options, &AA, perm_c, perm_r, &L, &U, &BB, &stat, &info);
		end2 = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds2 = end2-start2;
		std::cout << "Time taken to factorize using superLu with Reverese cuthill_reordering: " <<  (elapsed_seconds2).count() << std::endl;

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
		// Check for errors in factorization
		if (info == 0) {
		// Factorization successful
		// Timing for solving
			start = std::chrono::system_clock::now();
			dgstrs(NOTRANS, &L, &U, perm_c, perm_r, &BB, &stat, &info);
			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			//std::cout << "Time taken to solve Ax=B using SuperLU with Reverse cuthill re_ordering " << (elapsed_seconds).count() << std::endl;
		} else {
			std::cerr << "SuperLU factorization failed with error code: " << info << std::endl;
		}

		//cleaningup dynamically allocated variable.
		StatFree(&stat);
		Destroy_SuperMatrix_Store(&AA);
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
		Destroy_SuperMatrix_Store(&BB);
		A.clear();

		delete []perm_r;
		delete []perm_c;
	}
*/
/* ######################################################################################################
############################                BASIC SUPERLU           #####################################
###################################################################################################### */
	
	SuperMatrix AA, BB, L, U;
	
	dCreate_CompCol_Matrix(&AA, m, n, A.nnz, A.a.data(), A.lineindex.data(), A.columnptr.data(), SLU_NC, SLU_D, SLU_GE);
		
	dCreate_Dense_Matrix(&BB, m, nrhs, B.data(), m, SLU_DN, SLU_D, SLU_GE);
	
	superlu_options_t options;  

	// Set the default options:
	set_default_options(&options);	
	
	//options.ColPerm = METIS_AT_PLUS_A;  // graph partitioning and ordering tool that can be used to reorder the columns of the matrix.
	/* 
	options.IterRefine = SLU_DOUBLE;    // better choode double precision as the matrix is ill-conditioned
	options.Fact = DOFACT;              // instructs SuperLU to perform the mathematical operations necessary to break down the original matrix into LU
	options.Equil = YES;                // Equilibration scales the matrix to have unit row and column norms
	options.DiagPivotThresh = 1.0e-6;   // to tacke the ill-conditioned matrix
	options.SymmetricMode = YES;        // indicates that the input matrix is symmetric
	options.PivotGrowth = YES;          // enables tracking of pivot growth during the factorization
	options.ConditionNumber = NO;      // ndicates that SuperLU will estimate the condition number of the matrix. 
	options.PrintStat = YES;	    // print various statistics related to the factorization process and the solution. 
	*/
	

	SuperLUStat_t stat;             // this variable is used to store some statistics.
	StatInit(&stat);
	int info;                       // this variable contains informations on the succes or failure of the superLU fonction
	int *perm_c = new int [n];      // This variable stores the permutations applied on the columns of the matrix
	int *perm_r =  new int[m];      // This Variable stores the permutations applied on the row of the matrix


	 // Timing for factorization
	start = std::chrono::system_clock::now();
	dgssv(&options, &AA, perm_c, perm_r, &L, &U, &BB, &stat, &info);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Time taken to factorize suing SUPERLU: " <<  (elapsed_seconds).count() << std::endl;

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
	// Check for errors in factorization
	if (info == 0) {
	// Factorization successful
	// Timing for solving
		start = std::chrono::system_clock::now();
		dgstrs(NOTRANS, &L, &U, perm_c, perm_r, &BB, &stat, &info);
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "Time taken to solve Ax=B using SuperLU " << (elapsed_seconds).count() << std::endl;
	} else {
		std::cerr << "SuperLU factorization failed with error code: " << info << std::endl;
	}
	

	//cleaningup dynamically allocated variable.
	StatFree(&stat);
	Destroy_SuperMatrix_Store(&AA);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
	Destroy_SuperMatrix_Store(&BB);
	A.clear();

	delete []perm_r;
	delete []perm_c;
}
   



	
	
	
	

	
	
	
	

