#include <iostream>
#include <cmath>
#include <mpi.h>
#include <cassert>
#include <cstdlib>


int main(int argc, char** argv)
{
	int rank, size;
	int n = 240;
	
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
	double sum = 0, start_time, data_prep_duration, scatter_duration, calc_duration;
	
	if(rank == 0)
	{
		start_time = MPI_Wtime();
		double* data = new double[n];
		for(double i = 1; i <= n; i++)
			data[(int)i] = 1./pow(i, 2);
		data_prep_duration = MPI_Wtime() - start_time;
		MPI_Scatter(data, (int)n/size, MPI_DOUBLE, recv_data, n/size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		scatter_duration = MPI_Wtime() - start_time - data_prep_duration;
	} else {
		MPI_Scatter(NULL, (int)n/size, MPI_DOUBLE, recv_data, (int)n/size, MPI_DOUBLE, 0,  MPI_COMM_WORLD); 
	}
	
	for(int i = 0; i < n/size; i++)
		sum += recv_data[i];
	
	if(rank == 0)
	{
		calc_duration = MPI_Wtime() - start_time - data_prep_duration - scatter_duration;
		double pi;
		MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		pi = sqrt(6*pi);
		double duration = MPI_Wtime() - start_time;
		std::cout << pi << std::endl;
		std::cout << "Time: " << duration << std::endl;
		std::cout << "Data Prep: " << data_prep_duration << "\t Scatter: " << scatter_duration << "\t Calc: "<<  calc_duration
											<< "\tReduce: " << duration - data_prep_duration - scatter_duration - calc_duration  << std::endl;
		std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
	} else
	{
		MPI_Reduce(&sum, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	
	
	
	MPI_Finalize();
	return 0;
}

