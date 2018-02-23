#include <iostream>
#include <cmath>
#include <mpi.h>
#include <cassert>
#include <cstdlib>


#define MPITAG_GLOBAL_REDUCTION 100

double GlobalReductionSum(const double& res, const int rank, const int size)
{
	double s = res, buf=0;
	for(unsigned int d = 0; d < log2(size); d++)
	{
		int q =rank ^ (int)pow(2, d);
		//std::cout << rank << " sends " << s << " to " << q << std::endl;
		int err = MPI_Sendrecv(&s, 1, MPI_DOUBLE, q, MPITAG_GLOBAL_REDUCTION,
					&buf, 1, MPI_DOUBLE, q, MPITAG_GLOBAL_REDUCTION,
					MPI_COMM_WORLD, NULL);
		if(err != MPI_SUCCESS) std::cout << rank << ": Error " <<  err << std::endl;
		//std::cout << rank << " received " << buf << " from " << q << std::endl;
		s += buf;
	}
	return s;
}


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
			data[(int)i] = pow(-1, i-1)*(4.*pow(1./5., 2*i-1) - pow(1./239., 2*i-1))/(2*i-1);
		data_prep_duration = MPI_Wtime() - start_time;
		MPI_Scatter(data, (int)n/size, MPI_DOUBLE, recv_data, n/size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		scatter_duration = MPI_Wtime() - start_time - data_prep_duration;
	} else {
		MPI_Scatter(NULL, (int)n/size, MPI_DOUBLE, recv_data, (int)n/size, MPI_DOUBLE, 0,  MPI_COMM_WORLD); 
	}
	
	for(int i = 0; i < n/size; i++)
		sum += recv_data[i];
	
	if(rank == 0)
		calc_duration = MPI_Wtime() - start_time - data_prep_duration - scatter_duration;
	
	double pi;
	//MPI_Allreduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	pi = GlobalReductionSum(sum, rank, size);
	pi = 4*pi;
	
	if(rank == 0)
	{
		double duration = MPI_Wtime() - start_time;
		std::cout << pi << std::endl;
		std::cout << "Time: " << duration << std::endl;
		std::cout << "Data Prep: " << data_prep_duration << "\t Scatter: " << scatter_duration << "\t Calc: "<<  calc_duration
											<< "\tReduce: " << duration - data_prep_duration - scatter_duration - calc_duration  << std::endl;
		std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}

