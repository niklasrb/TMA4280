#include <cstdlib>
#include <mpi.h>
#include <cmath>
#include "../include/tma.h"

using namespace tma;


mesh<interval, 1> CreateIntervalMesh(real a, real b, int N) ;

int main(int argc, char**argv)
{
	
	int rank, size;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Op MPI_ERRSUM;
	MPI_Op_create( &tma::ErrorSum, true, &MPI_ERRSUM);
	
	uint N = 30;
	if(argc > 1) N = std::atoi(argv[1]); 

	distributedMesh<interval, 1> dm(CreateIntervalMesh(0., 1., N) , size);
	if(rank == 0) dm.dump();
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	// 1D Poisson Problem
	const double h = 1./N;
	DistributedVector x_i(dm, rank), b(dm, rank);
	DistributedSparseMatrix A(dm, N+1, rank);
	auto f = [](const double& x) { return std::sin(M_PI*x); };
	
	// Assemble
	
	for(uint i = A.RowRange().first; i < A.RowRange().second; i++) {
		A.set(i, i, 2.);
		if(i+1 < A.col())
			A.set(i, i+1, -1.);
		if(i > 0)
			A.set(i, i-1, -1.);
	}
	std::cout << A << std::endl;
	for(uint i = 0; i < dm.nverts(); i++) {
		if(dm.vertOwner(i) == (uint)rank) {
			b[i] = h*h * f(dm.geom()(i)[0] );
			x_i[i] = 0;
		}
	}
	tma::Sync(dm, x_i, rank);
	b.dump();
	x_i.dump();
	
	// Iteration
	DistributedVector x_i1(x_i);
	real globalNorm, localNorm, eps = 1e-6;
	int count = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	
	do {
		tma::JacobiStep(dm, rank, A, b, x_i);
		tma::Sync(dm, x_i, rank);
		localNorm =  (x_i + (-x_i1)).norm(); 
		MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_ERRSUM, MPI_COMM_WORLD);
		//if(rank == 0)
			//std::cout << count << ": " << globalNorm << " - " << x << std::endl;
		x_i1 = x_i; count++;
	}	while(globalNorm > eps && count < 1000);
	
	std::cout << "ps " << rank << ": " << x_i << std::endl;
	if(rank == 0)
		std::cout << "norm: " << globalNorm << " - count " << count << std::endl;
	
	DistributedVector x0(x_i);
	for(uint i = 0; i < dm.nverts(); i++) {
		if(dm.vertOwner(i) == (uint)rank) 
			x0[i] = f(dm.geom()(i)[0])/M_PI/M_PI;
	}
	localNorm = (x_i + (-x0)).norm();
	MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_ERRSUM, MPI_COMM_WORLD);
	auto res = tma::Collect(dm, x_i, rank, size);
	if(rank == 0)
		std::cout << res << std::endl << "expected " << x0 <<  std::endl << "residue " << (x_i + (-x0)) << " with norm " << globalNorm << std::endl;
	
	
	MPI_Finalize();
	return 0;
}

mesh<interval, 1> CreateIntervalMesh(real a, real b, int N) 
{
	mesh<interval, 1> m(N, N+1);
	for(uint i = 0; i < N; i++ ) {
		m.topo()(i)[0] = i;
		m.topo()(i)[1] = i+1;
	}
	for(uint i = 0; i < N+1; i++)
		m.geom()(i)[0] = a + (b-a)*double(i)/N;
	return m;
}


