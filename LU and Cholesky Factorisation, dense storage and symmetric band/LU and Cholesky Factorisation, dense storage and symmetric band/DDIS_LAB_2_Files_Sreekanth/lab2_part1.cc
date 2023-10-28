#include <algorithm>
#include <assert.h>
#include <iostream>
#include <fstream>

#include <random>
#include <chrono>
#include <functional>
#include <math.h>
#include "dmatrix_denseCM.h"
#include "dLU_denseCM.h"

int main(){

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.1, 1.);
	

	std::vector<int> input = {1000, 2000, 3000, 4000, 5000};

	std::ofstream outputFile("timing_results_Basic.txt");
    
	if (!outputFile.is_open()) {
		std::cerr << "Failed to open the file for writing." << std::endl;
		return 1;
	}

	for (int m: input) 
	{
	
		std::chrono::time_point<std::chrono::system_clock> start, end;
  		
  		start = std::chrono::system_clock::now();
		// This exemple generate m*m positive definite dense matrix A
		// Then factorise it, and solve the linear system Ax =B
		// Then we compute the residual AX-B to check if the factorisation is correct, 
		//  by computing norm2(Ax -B)
		//int m = 6;
		dmatrix_denseCM M(m,m);
		// generate random term in M. All terms are >=0.1 in M
		std::generate(M.begin(), M.end(), std::bind(dis,gen) );
		//A=MM^T -> This will insure that A is symetric positive definite : The factorisation LU always exist, as well as the cholesky  factorisation (A =LLT =UTU).

		dmatrix_denseCM A = M*transpose(M);

		// This build the factorisation of A, using the Basic algorithm.
		// IT will only work if you have properly implemented your factorisation.
		dLU_denseCM LU(A, factorpolicy (factorpolicy::Basic));
		// Below are call to factorisation using other version of the factorisation Algorithm.. You have to finish up Basic,L2,L3. Lapack is given.  
		//dLU_denseCM LU(A, factorpolicy(factorpolicy::L2));
		//dLU_denseCM LU(A, factorpolicy(factorpolicy::L3, 12));
		//dLU_denseCM LU(A, factorpolicy(factorpolicy::Lapack));


		// This will extract from the factorisation P, L and U, and build A2 = P^TLU, to check that your factor are corrects.
		// You don't want to do that with large matrix, this is just an illustration to help you check your result.

		dmatrix_denseCM L = LU.getL();
		dmatrix_denseCM U = LU.getU();
		dmatrix_denseCM P = LU.getP();
		/*
		std::cout << "A: " <<  A << std::endl;
		std::cout << "P: " << P << std::endl;
		std::cout << "L: " << L <<std::endl;
		std::cout << "U: " << U <<std::endl;
		*/
		//dmatrix_denseCM A2 = transpose(P)*L*U;
		
		//std::cout << "A2 = P^TLU : " << A2 <<std::endl;
		//std::cout << "error in the factorisation " <<  norm2(A-A2) << std::endl;
		
		dmatrix_denseCM B(m, 1);
		std::generate(B.begin(), B.end(), std::bind(dis,gen));
		dmatrix_denseCM X = LU.solve(B);
		
		end = std::chrono::system_clock::now();

		std::chrono::duration<double> duration = end - start;

		outputFile << "M: " << m << " Time: " << duration.count() << " seconds" << ", error: " << norm2(A*X - B) << std::endl;
		
		std::cout << "error in solving the system of size: " << m << " is: " <<  norm2(A*X - B) << " and time taken is: " << duration.count() << " seconds" << std::endl;
		
	}
		
	outputFile.close();
}
