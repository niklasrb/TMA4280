#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

double mach(unsigned int n)
{
	double sum = 0;
	for(unsigned int i = 1; i <= n; i++)
		sum += pow(-1, i-1)*(4.*pow(1./5., 2*i-1) - pow(1./239., 2*i-1))/(2*i-1);
	return 4.*sum;
}

bool unitTest(const double err = 1e-5)
{
	double res = mach(3);
	double exp = 4.*(4./5. - 1./239. - (4.*pow(1./5., 3) - pow(1./239., 3))/3. + (4.*pow(1./5., 5) - pow(1./239., 5))/5.) ; 
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
		ss << "n = 2^" << k << " = " << pow(2,k) << ": |pi - pi_n| = " << std::abs(M_PI - mach(pow(2,k))) << std::endl;
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
	} else {
		double pi = mach(n);
		std::cout << "n = " << n << "   pi ~= " << pi << std::endl << "Error: " << std::abs(M_PI - pi) << std::endl;
	}
	return 0; 
}
