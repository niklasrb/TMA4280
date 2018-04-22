#include <cstdlib>
#include <mpi.h>
#include <cmath>
#include "../include/tma.h"

using namespace tma;

#define MPI_SYNC 100


void Sync(DistributedVector& v, const distributedMesh<interval, 1>& dM, int rank);
void ErrorSum( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype);
template<class T, uint D>
void JacobiStep(const distributedMesh<T, D>& dm, uint rank, const DistributedSparseMatrix& A, const DistributedVector& b,  DistributedVector& x);

int main(int argc, char**argv)
{
	
	int rank, size;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Op MPI_ERRSUM;
	MPI_Op_create( &ErrorSum, true, &MPI_ERRSUM);
	
	uint N = 30;
	if(argc > 1) N = std::atoi(argv[1]); 
	
	mesh<interval, 1> m(N, N+1);
	for(uint i = 0; i < N; i++ ) {
		m.topo()(i)[0] = i;
		m.topo()(i)[1] = i+1;
	}
	for(uint i = 0; i < N+1; i++)
		m.geom()(i)[0] = double(i)/N;

	distributedMesh<interval, 1> dm(m, size);
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
	b.dump();
	x_i.dump();
	
	// Iteration
	DistributedVector x_i1(x_i);
	real globalNorm, localNorm, eps = 1e-8;
	int count = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	
	do {
		JacobiStep(dm, rank, A, b, x_i);
		Sync(x_i, dm, rank);
		localNorm =  (x_i + (-x_i1)).norm(); 
		MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_ERRSUM, MPI_COMM_WORLD);
		//if(rank == 0)
			//std::cout << count << ": " << globalNorm << " - " << x << std::endl;
		x_i1 = x_i; count++;
	}	while(globalNorm > eps && count < 100000);
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
	if(rank == 0)
		std::cout << "expected " << x0 <<  std::endl << "residue " << (x_i + (-x0)) << " with norm " << globalNorm << std::endl;
	
	
	MPI_Finalize();
	return 0;
}

template<class T, uint D>
void JacobiStep(const distributedMesh<T, D>& dm, uint rank, const DistributedSparseMatrix& A, const DistributedVector& b,  DistributedVector& x)
{	
	int v = 0;
	std::vector<real> s; s.resize(x.size());
	for(uint i = 0; i < dm.nverts(); i++)
	{
		if(dm.vertOwner(i) != rank)
			continue;
		
		s.at(v) = b[i];
		for(uint j = 0; j < A.col(); j++) {
			//std::cout << "A_" << i << "," << j << " = " << A(i, j) << std::endl;
			if(i!= j && A(i, j) != 0)
				s.at(v) -= A(i, j)*x.at(j);
		}
		s.at(v++) /= A(i, i);
	}
	v = 0;
	for(uint i = 0; i < dm.nverts(); i++) {
		if(dm.vertOwner(i) == rank)
			x[i] = s.at(v++);
	}	
}

void Sync(DistributedVector& v, const distributedMesh<interval, 1>& dM, int rank)
{
	double buf = 0;
	for(uint i = 0; i < dM.nverts(); i++) {
		if(dM.IsGhost(i)) {
			if(dM.vertOwner(i) == (uint)rank) {
				for(uint j = 0; j < dM.secondaryOwners(i).size(); j++)
					MPI_Send(&v[i], 1, MPI_DOUBLE, dM.secondaryOwners(i)[j], MPI_SYNC, MPI_COMM_WORLD);
			} else if( dM.secondaryOwners(i).contains(rank) ) {
				MPI_Recv(&buf, 1, MPI_DOUBLE, dM.vertOwner(i), MPI_SYNC, MPI_COMM_WORLD, NULL);
				v.update(i, buf);
			}
		}
	}
	
}

void ErrorSum( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
	for(uint i = 0; i < *len; i++) 
		((double*)inoutvec)[i] = std::sqrt(std::pow(((double*)invec)[i], 2) + std::pow(((double*)inoutvec)[i], 2));
}
