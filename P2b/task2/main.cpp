/* Niklas Becker
 * main.cpp 
 * task 2
 * 23.4.18
 */


#include "../include/tma.h"

using namespace tma;

mesh<triangle> CreateRectangularMesh(real x0, real y0, real Lx, real Ly, uint n);

int main(int argc, char**argv)
{
	
	int rank, size;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Op_create( &tma::ErrorSum, true, &MPI_ERRSUM);
	
	uint k = 3, N;
	if(argc > 1) k = std::atoi(argv[1]); 
	N = pow(2, k);
	
	// Create Mesh
	auto mesh = CreateRectangularMesh(0, 0, 1, 1, N);
	
	distributedMesh<triangle> dm(mesh, size);
	if(rank == 0)
		dm.dump();
	
	// Linearize Problem
	auto f = [] (const point<2>& x) { return (1.+2*M_PI*M_PI)*std::cos(M_PI*x(0))*std::sin(M_PI*x(1)); };
	
	DistributedSparseMatrix K(dm, dm.nverts(), rank), M(dm, dm.nverts(), rank), A(dm, dm.nverts(), rank);
	DistributedVector F(dm, rank);
	
	Assemble(dm, f, rank, K, M, F);
	A = K + M;
	Sync(dm, F, rank);
	if(rank == 0) {
		std::cout << "K = " << K << std::endl;
		std::cout << "M = " << M << std::endl;
		std::cout << "F = " << F << std::endl;
		std::cout << "A = " << A << std::endl;
	}
	
	// Solve
	DistributedVector u = ConjugateGradient(dm, A, F, rank, 1e-6, 10000);
	DistributedVector Au(dm, rank);
	std::cout << "u = " << u << std::endl;
	DistributedMatrixVectorProduct(dm, A, u, Au);
	real norm = (Au + (-F)).norm(); MPI_Allreduce(&norm, &norm, 1, MPI_DOUBLE, MPI_ERRSUM, MPI_COMM_WORLD);
	
	if(rank == 0)
		std::cout << "||Au - F|| = " << norm << std::endl;
	
	vector<real> u_coll = Collect(dm, u, rank, 0);
	
	
	// Compare
	if(rank == 0) {
		std::cout << "u_num = " << u_coll << std::endl;
		auto u_theoretical = [] (const point<2>& x) { return std::cos(M_PI*x(0))*std::sin(M_PI*x(1)); };
		auto u_num = [&u_coll, &dm] (const point<2>& x) 
		{
			triangle T; real s = 0;basefunctions<triangle> bf;
			for(uint c = 0; c < dm.ncells(); c++) {
				T = dm.physicalCell(c);
				for(uint i = 0; i < 3; i++) {
					s += u_coll.at(dm.topo()(c)(i))*bf.phi(i, T, x);
				}
			}
			return s;
		};
		auto L2integrand = [&u_theoretical, &u_num] (const point<2>& x) { return pow(u_theoretical(x) - u_num(x), 2); };
		real norm = 0;
		TWBQuadrature quad; 
		for(uint c = 0; c < dm.ncells(); c++) {
			norm += quad.Integral(dm.physicalCell(c), L2integrand);
		}
		norm = std::sqrt(norm);
		std::cout << "||u_theo(x) - u_num(x)||_2 = " << norm << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}

mesh<triangle> CreateRectangularMesh(real x0, real y0, real Lx, real Ly, uint n)
{
	mesh<triangle> m( 2*(n-1)*(n-1), n*n);
	for(uint k = 0; k < n*n; k++)
	{
		uint i = k % n; uint j = (k - i)/n;
		m.geom()(k) = point<2>({x0 + i*(Lx/real(n-1)), y0 + j*(Ly/real(n-1))});
		if(i < n-1 && j < n-1) {
			m.topo()(j*2*(n-1) + 2*i) = cell<triangle::N()>({j*n + i, (j+1)*n + i, j*n + i + 1}); 
			m.topo()(j*2*(n-1) + 2*i + 1) = cell<triangle::N()>({(j+1)*n + i+1, (j+1)*n + i, j*n + i + 1}) ; 
		}
	}
	return m;
}
