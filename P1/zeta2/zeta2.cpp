#include <iostream>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <ctime>


inline double RZ(double i)
{
	return 1./pow(i, 2);
}

int main(int argc, char** argv)
{
	int n = 3;
	if(argc > 1) n = atoi(argv[1]);
	
	double pi = 0, duration = clock();
	
	
	#pragma omp parallel for reduction(+:pi)
	for(int i = 1; i <= n; i++)  {
		pi += RZ((double)i);
	}
	pi = sqrt(6*pi);
	duration = (clock() - duration)/CLOCKS_PER_SEC;
	
	std::cout << pi << std::endl;
	std::cout << duration << std::endl;
	std::cout << std::abs(M_PI - pi) << std::endl;
	
	return 0;
}
