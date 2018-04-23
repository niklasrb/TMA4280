#ifndef TMA_JACOBI_H_
#define TMA_JACOBI_H_

#include <tma/types.h>
#include <cmath>
#include <mpi.h>

namespace tma
{
#define MPI_SYNC 100
	
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
			if(i!= j && A(i, j) != 0)
				s.at(v) -= A(i, j)*x.at(j);
		}
		s.at(v) /= A(i, i); v++;
	}
	v = 0;
	for(uint i = 0; i < dm.nverts(); i++) {
		if(dm.vertOwner(i) == rank) {
			x[i] = s.at(v); v++;
		}
	}	
}

template<class T, uint D>
void Sync(const distributedMesh<T, D>& dM, DistributedVector& v, int rank)
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
template<class T, uint D>
vector<real> Collect(const distributedMesh<T, D>& dM, const DistributedVector& v, int rank, int size)
{
	double buf = 0;
	vector<real> cv;
	for(uint i = 0; i < dM.nverts(); i++) {
		if(dM.vertOwner(i) == rank) 
			buf = v[i];
		MPI_Bcast(&buf, 1, MPI_DOUBLE,  dM.vertOwner(i), MPI_COMM_WORLD);
		cv.push_back(buf);
	}
	return cv;
}

void ErrorSum( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
	for(uint i = 0; i < *len; i++) 
		((double*)inoutvec)[i] = std::sqrt(std::pow(((double*)invec)[i], 2) + std::pow(((double*)inoutvec)[i], 2));
}
} /* namespace tma */

#endif /* TMA_JACOBI_H_ */
