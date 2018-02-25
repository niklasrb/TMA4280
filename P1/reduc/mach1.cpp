#include <iostream>
#include <cmath>
#include <mpi.h>
#include <cassert>
#include <cstdlib>


#define MPITAG_GLOBAL_REDUCTION 100

double GlobalReductionSum(double value, const int rank, const int size)
{
	double buf=0;
	int q, err;
	for(unsigned int d = 0; d < log2(size); d++)
	{
		q = rank ^ (int)pow(2, d);
		//std::cout << rank << " sends " << s << " to " << q << std::endl;
		err = MPI_Sendrecv(&value, 1, MPI_DOUBLE, q, MPITAG_GLOBAL_REDUCTION,
					&buf, 1, MPI_DOUBLE, q, MPITAG_GLOBAL_REDUCTION,
					MPI_COMM_WORLD, NULL);
		if(err != MPI_SUCCESS) 
			std::cout << rank << ": Error " <<  err << std::endl;
		else
			value += buf;
		//std::cout << rank << " received " << buf << " from " << q << std::endl;
	}
	return value;
}

inline double mach(double i)
{
	return ((int)i%2==0 ? -1. : 1.)*(4.*pow(1./5., 2*i-1) - pow(1./239., 2*i-1))/(2.*i-1.);
}

int main(int argc, char** argv)
{
	int rank, size;
	unsigned int n = 1, reduction_mode = 0;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(argc > 1) n = atoi(argv[1]);
	if(argc > 2) reduction_mode = atoi(argv[2]);
	
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&reduction_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(n % size != 0 || floor(log2(size)) != log2(size)) {
		if(rank == 0) std::cout << "Argument Error" << std::endl;
		MPI_Finalize();
		return 0;
	}
	
	double sum = 0;
	double start_time, calc_duration, reduce_duration;
	unsigned int left = (unsigned)rank*n/(unsigned)size +1;
	unsigned int right = ((unsigned)rank+1)*n/(unsigned)size +1;
	
	if(rank == 0)
		start_time = MPI_Wtime();
	
	for(unsigned int i = left; i <= right; i++)
		sum += mach(i);
	
	if(rank == 0)
		calc_duration = MPI_Wtime() - start_time;
	MPI_Barrier(MPI_COMM_WORLD);
	
	reduce_duration = MPI_Wtime();
	double pi;
	switch(reduction_mode) {
		case 0: MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			break;
		case 1: MPI_Allreduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
			break;
		case 2: pi =  GlobalReductionSum(sum, rank, size);
			break;
	}
	pi = 4*pi;
	if(rank == 0) {
		double duration = MPI_Wtime() - start_time;
		reduce_duration = MPI_Wtime() - reduce_duration;
		std::cout << pi << std::endl;
		std::cout << "Time: " << duration << "  - Calc: "<<  calc_duration	<< "\tReduce: " << duration - calc_duration  << std::endl;
		std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
		std::cout << "Reduction mode " << reduction_mode << std::endl;
	}
	
	
	
	MPI_Finalize();
	return 0;
}
