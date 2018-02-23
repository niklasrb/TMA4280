#include <iostream>
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <cstdlib>
#include <chrono>


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
	
	double* recv_data = new double[n/size];
	double sum = 0;
	std::chrono::time_point<std::chrono::high_resolution_clock> start, data_prep, scatter, calc, end;
	
	if(rank == 0)
	{
		start = std::chrono::high_resolution_clock::now();
		double* data = new double[n];
		for(double i = 1; i <= n; i++)
			data[(int)i] = pow(-1, i-1)*(4.*pow(1./5., 2*i-1) - pow(1./239., 2*i-1))/(2*i-1);
		data_prep = std::chrono::high_resolution_clock::now();
		MPI_Scatter(data, (int)n/size, MPI_DOUBLE, recv_data, n/size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		scatter = std::chrono::high_resolution_clock::now();
	} else {
		MPI_Scatter(NULL, (int)n/size, MPI_DOUBLE, recv_data, (int)n/size, MPI_DOUBLE, 0,  MPI_COMM_WORLD); 
	}
	
	#pragma omp parallel for reduction(+:sum)
	for(int i = 0; i < n/size; i++)
		sum += recv_data[i];
	
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
		std::cout << "Data Prep: " << (std::chrono::duration<double>(data_prep - start)).count() << " s" << "\t Scatter: " << (std::chrono::duration<double>(scatter - data_prep)).count() << " s" 
											<< "\t Calc: "<< (std::chrono::duration<double>(calc - scatter)).count() << " s"
											<< "\tReduce: " << (std::chrono::duration<double>(end - calc)).count() << " s"  << std::endl;
		std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}

