#ifndef TMA_PARALLEL_H
#define TMA_PARALLEL_H

#include <tma/types.h>
#include <mpi.h>

namespace tma
{
	
MPI_Op MPI_ERRSUM;
#define MPI_SYNC 100

void ErrorSum( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
	for(uint i = 0; i < *len; i++) 
		((double*)inoutvec)[i] = std::sqrt(std::pow(((double*)invec)[i], 2) + std::pow(((double*)inoutvec)[i], 2));
}

template<class T>
void Sync( const distributedMesh<T>& dM, DistributedVector& v, uint rank)
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
template<class T>
DistributedVector& DistributedMatrixVectorProduct(const distributedMesh<T>& dM, const DistributedSparseMatrix& A, const DistributedVector& x, DistributedVector& res)
{
	real s;
	for(uint i = A.RowRange().first; i < A.RowRange().second; i++) {
		//if(dM.vertOwner(i) != rank) continue;
		s = 0;
		for(uint j = 0; j < A.col(); j++) {
			if(A(i, j) != 0)
				s += A(i, j)*x[j];
		}
		res[i] = s;
	}
	return res;
}

template<class T>
DistributedVector ConjugateGradient(const distributedMesh<T>& dm, const DistributedSparseMatrix& A, const DistributedVector& b, uint rank, real eps = 1e-5, uint maxIt = 1)
{
	DistributedVector x, r, rprev, e, Ae, buf;
	real alpha, norm, eAe; uint i = 0 ;
	Sync(dm, x, rank); 
	r = b - DistributedMatrixVectorProduct(dm, A, x, buf); Sync(dm, r, rank);
	std::cout << "r = " << r << std::endl;
	e = r; DistributedMatrixVectorProduct(dm, A, x, Ae); Sync(dm, Ae, rank);
	std::cout << "A*e = " << Ae << std::endl;
	eAe = e*Ae; MPI_Allreduce(&eAe, &eAe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	alpha = r*b/eAe; MPI_Allreduce(&alpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	std::cout << "alpha " << alpha << std::endl;
	x = alpha*e; Sync(dm, x, rank);
	std::cout << "x = " << x << std::endl;
	do {
		rprev = r;
		r = rprev - alpha*Ae;
		Sync(dm, r, rank);
		e = r - (1./eAe)*e*DistributedMatrixVectorProduct(dm, A, rprev, buf);
		Sync(dm, e, rank);
		DistributedMatrixVectorProduct(dm, A, e, Ae);
		Sync(dm, Ae, rank);
		eAe = e*Ae;
		MPI_Allreduce(&eAe, &eAe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		alpha = e*b/eAe;
		MPI_Allreduce(&alpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		x = x + alpha*e;
		Sync(dm, x, rank);
		norm = alpha*e.norm();
		MPI_Allreduce(&norm, &norm, 1, MPI_DOUBLE, MPI_ERRSUM, MPI_COMM_WORLD);
		
	} while( norm > eps && ++i < maxIt );
	return x;
}
	


} // namespace tma


#endif // TMA_PARALLEL_H
