#ifndef TMA_PARALLEL_H
#define TMA_PARALLEL_H

#include <tma/types.h>
#include <mpi.h>

namespace tma
{
	
MPI_Op MPI_ERRSUM;
#define MPI_SYNC 100
#define MPI_COLLECT 101

void ErrorSum( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
	for(int i = 0; i < *len; i++) 
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
vector<real> Collect(const distributedMesh<T> dm, const DistributedVector& a, uint rank, uint root)
{
	vector<real> b; b.resize(dm.nverts()); real buf;
	for(uint v = 0; v < dm.nverts(); v++) {
		if(dm.vertOwner(v) == rank) {
			if(rank == root)  
				b[v] = a[v];
			else {
				buf = a[v];
				MPI_Send(&buf, 1, MPI_DOUBLE, root, MPI_COLLECT, MPI_COMM_WORLD);
			} 
		} else if(root == rank) {
			MPI_Recv(&b[v], 1, MPI_DOUBLE, dm.vertOwner(v), MPI_COLLECT, MPI_COMM_WORLD, NULL);
		}
	}
	return b;
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
DistributedVector ConjugateGradient(const distributedMesh<T>& dm, const DistributedSparseMatrix& A, const DistributedVector& b, uint rank, real eps = 1e-5, uint maxIt = 100)
{
	DistributedVector x(dm, rank), r(dm, rank), rprev(dm, rank), e(dm, rank), Ae(dm, rank), buf(dm, rank);
	real alpha, beta, norm, eAe, rr; uint i = 0 ;
	Sync(dm, x, rank); 
	r = b - DistributedMatrixVectorProduct(dm, A, x, buf); Sync(dm, r, rank);
	//std::cout << rank << ": r = " << r << std::endl;
	e = r; 
	do {
		DistributedMatrixVectorProduct(dm, A, e, Ae); Sync(dm, Ae, rank);
		
		eAe = e*Ae; MPI_Allreduce(&eAe, &eAe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		alpha = r*r/eAe; MPI_Allreduce(&alpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		x = x + alpha*e; Sync(dm, x, rank);
		//std::cout << rank << ": x = " << x << std::endl;
		
		rprev = r;
		r = r - alpha*Ae; Sync(dm, r, rank);
		//std::cout << rank << ": r = " << r << std::endl;
		
		norm = r.norm(); MPI_Allreduce(&norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if(norm < eps) continue;
		
		rr = rprev*rprev; MPI_Allreduce(&rr, &rr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		beta = (1/rr)*r*r; MPI_Allreduce(&beta, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		e = r + beta*e; Sync(dm, e, rank);		
	} while( norm > eps && ++i < maxIt );
	return x;
}
	


} // namespace tma


#endif // TMA_PARALLEL_H
