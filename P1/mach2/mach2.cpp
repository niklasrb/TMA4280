#include <iostream>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <chrono>


inline double mach(double i)
{
	return pow(-1, i-1)*(4.*pow(1./5., 2*i-1) - pow(1./239., 2*i-1))/(2*i-1);
}

int main(int argc, char** argv)
{
	int n = 3;
	if(argc > 1) n = atoi(argv[1]);
	
	double pi = 0;
	auto start = std::chrono::high_resolution_clock::now();
	
	
	#pragma omp parallel for reduction(+:pi)
	for(int i = 1; i <= n; i++)  {
		pi += mach((double)i);
	}
	pi = 4*pi;
	
	auto end = std::chrono::high_resolution_clock::now();
	
	std::cout << " pi = " << pi << "\t n = " << n  << std::endl;
	std::cout << "Time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;
	std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
	
	return 0;
}
