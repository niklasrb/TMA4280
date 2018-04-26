#include <cstdlib>
#include <mpi.h>
#include <cmath>
#include "../include/tma.h"

using namespace tma;


mesh<interval> CreateIntervalMesh(real a, real b, uint N) ;

int main(int argc, char**argv)
{
	
	int rank, size;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Op MPI_ERRSUM;
	MPI_Op_create( &tma::ErrorSum, true, &MPI_ERRSUM);
	
	uint N = 8;
	if(argc > 1) N = std::atoi(argv[1]); 

	distributedMesh<interval> dm(CreateIntervalMesh(0., 1., N) , size);
	if(rank == 0) dm.dump();
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	// 1D Poisson Problem
	DistributedVector u_exp(dm, rank), f_exp(dm, rank);
	DistributedSparseMatrix L(dm, dm.nverts(), rank), M(dm, dm.nverts(), rank);
	auto f = [](const double& x) { return std::sin(M_PI*x); };
	
	
	
	// Assemble
	

	
	MPI_Finalize();
	return 0;
}

mesh<interval> CreateIntervalMesh(real a, real b, uint N) 
{
	mesh<interval> m(N, N+1);
	for(uint i = 0; i < N; i++ ) {
		m.topo()(i)(0) = i;
		m.topo()(i)(1) = i+1;
	}
	for(uint i = 0; i < N+1; i++)
		m.geom()(i)(0) = a + (b-a)*double(i)/N;
	return m;
}


