#include <iostream>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <chrono>


inline double RZ(double i)
{
	return 1./pow(i, 2);
}

int main(int argc, char** argv)
{
	int n = 3;
	if(argc > 1) n = atoi(argv[1]);
	
	double pi = 0;
	auto start = std::chrono::high_resolution_clock::now();
	
	
	#pragma omp parallel for reduction(+:pi)
	for(int i = 1; i <= n; i++)  {
		pi += RZ((double)i);
	}
	pi = sqrt(6*pi);
	
	auto end = std::chrono::high_resolution_clock::now();
	
	std::cout << "pi = " << pi << "\t n = " << n << std::endl;
	std::cout << "Time: " << std::chrono::duration<double, std::milli>(end-start).count() << " ms" << std::endl;
	std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
	
	return 0;
}
