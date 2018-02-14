#include <iostream>
#include <cmath>
#include <cstdlib>

double RZ(unsigned int n)
{
	double sum = 0;
	for(unsigned int i = 1; i <= n; i++)
		sum += 1./(i*i);
	return sum;
}


int main(int argc, char** argv)
{
	unsigned int n = 1;
	if(argc > 1) n = atoi(argv[1]);
	double res = RZ(n);
	std::cout << "n = " << n << "  RZ(n) = " << res << std::endl << "pi ~= " << sqrt(6*res) << std::endl;
	return 0; 
}
