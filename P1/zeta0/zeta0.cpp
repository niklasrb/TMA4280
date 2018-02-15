#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

long double RZ(unsigned int n)
{
	long double sum = 0;
	for(double  i = 1; i <= n; i++)
		sum += 1./pow(i, 2);
	return sqrt(6.*sum);
}

bool unitTest(const double err = 1e-4)
{
	double res = RZ(3);
	double exp = sqrt(6.*(1. + 1./4. + 1./9.));
	bool compare = std::abs(res - exp) < err;
	std::cout << "for n=3: " << res << std::endl;
	std::cout << "compared to " << exp << std::endl;
	std::cout << "Unit test completed " << (compare ? "succesfully": "unsucessfully") << std::endl;
	return compare;
}

void verificationTest(std::fstream& f, unsigned int maxK = 24)
{
	std::stringstream ss;
	for(unsigned int k = 1; k <= maxK; k++) {
		ss << "n = 2^" << k << " = " << pow(2,k) << ": |pi - pi_n| = " << std::abs(M_PI - RZ(pow(2,k))) << std::endl;
		std::cout << ss.str();
		if(f.good()) f << ss.str();
		ss.str(std::string());
	}
}

int main(int argc, char** argv)
{
	int n = 1;
	if(argc > 1) n = atoi(argv[1]);
	if(n == 0) 
		unitTest();
	else if(n == -1) {
		std::fstream f("verificationTest.txt", std::fstream::out | std::fstream::trunc);
		verificationTest(f);
		f.close();
	} else 
		std::cout << "n = " << n << "   pi ~= " << RZ(n) << std::endl;
	return 0; 
}
