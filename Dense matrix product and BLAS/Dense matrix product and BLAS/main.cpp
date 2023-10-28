// This file could be the starting point for your own program

#include <iostream>
#include <random>
#include <algorithm>
#include<functional>
#include <chrono>
#include <fstream>
#include "dmatrix_denseCM.h" 
#include <fstream>

using mat = dmatrix_denseCM;
using time_point = std::chrono::time_point<std::chrono::system_clock>;
using sclock =  std::chrono::system_clock;
using time_s = std::chrono::duration<double> ;


int main()
{

	int a[] = {10, 100, 1000, 2000, 3000, 4000};
	int size = sizeof(a)/sizeof(a[0]);
	
	std::fstream output_file;  //file handle
	
	output_file.open("output.txt", std::ios::out);
	
	for(int i = 0; i < size; ++i)
	{
		bool verbose = false;
		// select the size of your matrix below
		const size_t m = a[i];
		const size_t k = a[i];
		const size_t n = a[i];
		// Build 2 square matrix of size nxn
		mat A(m,k);
		mat B(k,n);

		// Setting up the random number generator
		std::random_device rd;
		std::mt19937 gen(rd());  // this is the generator
		std::uniform_real_distribution<> dis(0., 1.); // this 
		// Fill up the Matrix D with random number

		for (auto &a : A) a = dis(gen);
		for (auto &b : B) b = dis(gen);

		if (verbose) std::cout << A << std::endl;
		if (verbose) std::cout << B << std::endl;
		// compute A*B, measure the time it take and compute the computation speed.
					
		time_point tstart1 = sclock::now();
		mat C1 = mulV1(A,B);
		time_point tend1 = sclock::now();
		time_s time_mult1 = tend1 - tstart1;
		std::cout << m << "x" << n << " prod computed in time_mult " << time_mult1.count()<< "s" << std::endl;
		double r1 = 2*m*k*n/time_mult1.count()/1.e+9 ;
		std::cout << " Gflop/s "<< r1 << std::endl;
		if (verbose) std::cout << C1 << std::endl;	
		
		if(output_file.is_open())
		{
			output_file << "mulV1 " << a[i] << " " << time_mult1.count()  << "," << r1 << std::endl;
		}else
		{
			throw -1;
		}
		
		time_point tstart2 = sclock::now();
		mat C2 = mulV2(A,B);
		time_point tend2 = sclock::now();
		time_s time_mult2 = tend2 - tstart2;
		std::cout << m << "x" << n << " prod computed in time_mult " << time_mult2.count()<< "s" << std::endl;
		double r2 = 2*m*k*n/time_mult2.count()/1.e+9;
		std::cout << " Gflop/s "<< r2 << std::endl;
		if (verbose) std::cout << C2 << std::endl;
		
		if(output_file.is_open())
		{
			output_file << "mulV2 " << a[i] << " " << time_mult2.count() << "," << r2 << std::endl;
		}else
		{
			throw -1;
		}
			
		time_point tstart3 = sclock::now();
		mat C3 = muldgemm(A,B);
		time_point tend3 = sclock::now();
		time_s time_mult3 = tend3 - tstart3;
		std::cout << m << "x" << n << " prod computed in time_mult " << time_mult3.count()<< "s" << std::endl;
		double r3 =  2*m*k*n/time_mult3.count()/1.e+9;
		std::cout << " Gflop/s "<< r3 << std::endl;
		if (verbose) std::cout << C3 << std::endl;
		
		if(output_file.is_open())
		{
			output_file << "dgemm " << a[i] << " " << time_mult3.count() << "," << r3 << std::endl;
		}else
		{
			throw -1;
		}
	}
}

