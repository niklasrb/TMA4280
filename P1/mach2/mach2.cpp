#include <iostream>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <ctime>


inline double mach(double i)
{
	return pow(-1, i-1)*(4.*pow(1./5., 2*i-1) - pow(1./239., 2*i-1))/(2*i-1);
}

int main(int argc, char** argv)
{
	int n = 3;
	if(argc > 1) n = atoi(argv[1]);
	
	double pi = 0, duration = clock();
	
	
	#pragma omp parallel for reduction(+:pi)
	for(int i = 1; i <= n; i++)  {
		pi += mach((double)i);
	}
	pi = 4*pi;
	duration = (clock() - duration)/CLOCKS_PER_SEC;
	
	std::cout << pi << std::endl;
	std::cout << "Time: " << duration << std::endl;
	std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
	
	return 0;
}
