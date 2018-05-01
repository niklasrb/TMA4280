
#include "../include/tma.h"

using namespace tma;

void AssembleHelmholtz(const distributedMesh<interval>& dm, real (*f)(const point<1>&), uint rank, DistributedSparseMatrix& StiffnessMatrix, DistributedSparseMatrix& MassMatrix, DistributedVector& ForceVector);
void ExpandFunction(const distributedMesh<interval>& dm, const std::function<real(const point<1>&)>& f, uint rank, DistributedVector& f_expanded);
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
	//if(rank == 0) 
	//	dm.dump();
	if(rank == 0)
		std::cout << "N = " << N << " size = " << size << std::endl;
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	// 1D Poisson Problem
	DistributedVector u_exp(dm, rank), f_exp(dm, rank), res(dm, rank);
	DistributedSparseMatrix L(dm, dm.nverts(), rank), M(dm, dm.nverts(), rank);
	auto f = [](const point<1>& x) { return std::sin(M_PI*x(0)); };
	
	// Assemble
	AssembleHelmholtz(dm, f, rank, L, M, f_exp);
	
	Sync(dm, f_exp, rank);
	
	/*if(rank == 0) {
		std::cout << "L = " << L << std::endl;
		std::cout << "M = " << M << std::endl;
		std::cout << "f = " << f_exp << std::endl;
	}*/
	// Test
	auto u = [](const point<1>& x) { return std::sin(M_PI*x(0))/M_PI/M_PI; };
	ExpandFunction(dm, u, rank, u_exp);
	Sync(dm, u_exp, rank);
	
	//if(rank == 0)
	//	std::cout << "u = " << u_exp << std::endl;
	
	DistributedMatrixVectorProduct(dm, M+L, u_exp, res);
	if(rank == 0)
		std::cout << res << std::endl;

	real norm = (res - u_exp).norm(); MPI_Allreduce(&norm, &norm, 1, MPI_DOUBLE, MPI_ERRSUM, MPI_COMM_WORLD);
	if(rank == 0)
		std::cout << "||(M+L)u - f||_2 = " << norm << std::endl;
	
	MPI_Finalize();
	return 0;
}

void AssembleHelmholtz(const distributedMesh<interval>& dm, real (*f)(const point<1>&), uint rank, DistributedSparseMatrix& StiffnessMatrix, DistributedSparseMatrix& MassMatrix, DistributedVector& ForceVector)
{
	interval I; real Isize;
	uint globJ;
	basefunctions<interval> bf;
	GaussianQuadrature quad;
	
	#ifdef OPENMP
	#pragma omp parallel for private(I, Isize, globJ, bf, quad) schedule(guided)
	#endif
	for(uint v = 0; v < dm.nverts(); v++) {		// go through all vertices
		if(dm.vertOwner(v) != rank) continue;
		if(dm.geom()(v)(0) == 0. || dm.geom()(v)(0) == 1.) {	// boundary conditions
			StiffnessMatrix.set(v, v, 1.);
			ForceVector[v] = 0;
			continue;
		}
		for(uint c = 0; c < dm.ncells(); c++) {	// find cells the vertex belongs to
			if(dm.topo()(c).contains(v)) {
				I = dm.physicalCell(c);
				Isize = I.area();
				uint i = dm.topo()(c).find(v);
				for(uint j = 0; j < 2; j++) {	// find other vertices in cell
					globJ = dm.topo()(c)(j);
					StiffnessMatrix.set(v, globJ, StiffnessMatrix(v, globJ) + pow(-1., i)*pow(-1., j)/Isize/1.  );
					MassMatrix.set(v, globJ, MassMatrix(v, globJ) + quad.Integral(I, [&bf, i, j, &I] (const point<1>& x) { return bf.phi(i, I, x)*bf.phi(j, I, x); }));
				}
				ForceVector[v] += quad.Integral(I, [&bf, &f, i, &I] (const point<1>& x) { return f(x)*bf.phi(i, I, x);  });
			}
		}	
		
	}
	
}

void ExpandFunction(const distributedMesh<interval>& dm, const std::function<real(const point<1>&)>& f, uint rank, DistributedVector& f_expanded)
{
	interval I;
	basefunctions<interval> bf;
	GaussianQuadrature quad;
	
	for(uint v = 0; v < dm.nverts(); v++) {		// go through all vertices
		if(dm.vertOwner(v) != rank) continue;
		f_expanded[v] = 0.;
		for(uint c = 0; c < dm.ncells(); c++) {	// find cells the vertex belongs to
			if(dm.topo()(c).contains(v)) {
				I = dm.physicalCell(c);
				uint i = dm.topo()(c).find(v);
				f_expanded[v] += quad.Integral(I, [&bf, &f, i, &I] (const point<1>& x) { return f(x)*bf.phi(i, I, x);  });
			}
		}	
	}
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


