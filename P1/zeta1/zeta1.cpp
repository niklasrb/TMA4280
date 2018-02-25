#include <iostream>
#include <cmath>
#include <mpi.h>
#include <cassert>
#include <cstdlib>


inline double RZ(double i)
{
	return 1./pow(i, 2);
}

int main(int argc, char** argv)
{
	int rank, size;
	unsigned int n = 1;
	
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
	double start_time, calc_duration;
	unsigned int left = (unsigned)rank*n/(unsigned)size +1;
	unsigned int right = ((unsigned)rank+1)*n/(unsigned)size +1;
	
	if(rank == 0)
		start_time = MPI_Wtime();
	
	for(unsigned int i = left; i <= right; i++)
		sum += RZ(i);
	
	//std::cout << rank << ": from " << left << " till " << right << "  - " << sum << std::endl;
	
	if(rank == 0)
	{
		calc_duration = MPI_Wtime() - start_time;
		double pi;
		MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		pi = sqrt(6*pi);
		double duration = MPI_Wtime() - start_time;
		std::cout << pi << std::endl;
		std::cout << "Time: " << duration << "  - Calc: "<<  calc_duration	<< "\tReduce: " << duration - calc_duration  << std::endl;
		std::cout << "Error: " << std::abs(M_PI - pi) << std::endl;
	} else
	{
		MPI_Reduce(&sum, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	
	
	
	MPI_Finalize();
	return 0;
}

