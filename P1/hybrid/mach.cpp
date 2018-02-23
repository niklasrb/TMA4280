#include <iostream>
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <cstdlib>
#include <chrono>

inline double mach(double i)
{
	return pow(-1, i-1)*(4.*pow(1./5., 2*i-1) - pow(1./239., 2*i-1))/(2*i-1);
}


int main(int argc, char** argv)
{
	int rank, size;
	int n = 1;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(argc > 1) n = atoi(argv[1]);
	
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(n % size != 0 || floor(log2(size)) != log2(size)) {
		if(rank == 0) std::cout << "Argument Error" << std::endl;
		MPI_Finalize();
		return 0;
	}
	
	double sum = 0;
	std::chrono::time_point<std::chrono::high_resolution_clock> start, calc, end;
	
	if(rank == 0)
		start = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel for reduction(+:sum)
	for(int i = rank*n/size + 1; i <= (rank+1)*n/size; i++)
		sum += mach(i);
	
	if(rank == 0)
		calc = std::chrono::high_resolution_clock::now();
	
	double pi;
	MPI_Allreduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	pi = 4*pi;
	
	if(rank == 0)
	{
		end = std::chrono::high_resolution_clock::now();
		std::cout << pi << std::endl;
		std::cout << "Time: " << (std::chrono::duration<double>(end - start)).count() << " s" << std::endl;
		std::cout << "\t Calc: "<< (std::chrono::duration<double>(calc - start)).count() << " s"
											<< "\tReduce: " << (std::chrono::duration<double>(end - calc)).count() << " s"  << std::endl;
		std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}

