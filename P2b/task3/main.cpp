/* Niklas Becker
 * main.cpp 
 * task 3
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
	
	uint  N = 5, steps = 1;
	if(argc > 1) N = std::atoi(argv[1]); 
	if(argc > 2) steps = std::atoi(argv[2]); 
	
	// Create Mesh
	auto mesh = CreateRectangularMesh(0, 0, 1, 1, N);
	
	distributedMesh<triangle> dm(mesh, size);
	//if(rank == 0)
	//	dm.dump();
	
	// Linearize Problem
	real nu = 1.;
	real h = 0.01;
	auto u_tilde = [nu] (const point<2>& x, const real& t) { return exp(-nu*nu*t) * sin(sqrt(nu)*x(0) + sqrt(nu)*x(1)); };
	
	DistributedSparseMatrix K(dm, dm.nverts(), rank), M(dm, dm.nverts(), rank), S(dm, dm.nverts(), rank);
	
	MPI_Barrier(MPI_COMM_WORLD);
	auto startassembly = std::chrono::steady_clock::now();
	AssembleUnsteadyHeatProblem(dm, rank, K, M, S);
	MPI_Barrier(MPI_COMM_WORLD);
	auto endassembly= std::chrono::steady_clock::now();

	/*{
		std::cout << "K = " << K << std::endl;
		std::cout << "M = " << M << std::endl;
		std::cout << "S = " << S << std::endl;
	}*/
	
	if(rank == 0) {
		std::cout << "Unsteady Heat Equation" << std::endl;
		std::cout << "N = " << N << " ncells() = " << dm.ncells() << std::endl;
		std::cout << "size = " << size << std::endl;
		#ifdef OPENMP
		std::cout << "omp is used with " << omp_get_max_threads() << " threads" << std::endl;
		#else
		std::cout << "omp is not used" << std::endl;
		#endif
		std::cout << "Assembly took " << std::chrono::duration <double, std::milli> (endassembly - startassembly).count() << "ms" << std::endl;
	}
	
	// Solve
	MPI_Barrier(MPI_COMM_WORLD);
	auto startbe = std::chrono::steady_clock::now();
	DistributedVector u = BackwardsEuler(dm, M, K, S, u_tilde, nu, h, 0., rank, steps);
	MPI_Barrier(MPI_COMM_WORLD);
	auto endbe= std::chrono::steady_clock::now();
	
	//std::cout << "u = " << u << std::endl;
	
	vector<real> u_coll = Collect(dm, u, rank, 0);
	
	
	// Compare
	if(rank == 0) {
		//std::cout << "u_num = " << u_coll << std::endl;
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
		real t = h*steps;
		auto L2integrand = [&u_tilde, &u_num, &t] (const point<2>& x) { return pow(u_tilde(x, t) - u_num(x), 2); };
		real norm = 0;
		TWBQuadrature quad; 
		for(uint c = 0; c < dm.ncells(); c++) {
			norm += quad.Integral(dm.physicalCell(c), L2integrand);
		}
		norm = std::sqrt(norm);
		std::cout << "solving took " << std::chrono::duration <double, std::milli> (endbe - startbe).count() << "ms for " << steps << " steps" << std::endl;
		std::cout << "after t = " << t << " the constructed function has ||u_num(x) - u_tilde(x)||_2 = " << norm << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}

mesh<triangle> CreateRectangularMesh(real x0, real y0, real Lx, real Ly, uint n)
{
	mesh<triangle> m( 2*(n-1)*(n-1), n*n);
	for(uint k = 0; k < n*n; k++) // global index
	{
		uint i = k % n; uint j = (k - i)/n;	// get (i,j) index
		m.geom()(k) = point<2>({x0 + i*(Lx/real(n-1)), y0 + j*(Ly/real(n-1))});
		if(i < n-1 && j < n-1) {	// create triangles to the upper right of this vertex
			m.topo()(j*2*(n-1) + 2*i) = cell<triangle::N()>({j*n + i, (j+1)*n + i, j*n + i + 1}); 
			m.topo()(j*2*(n-1) + 2*i + 1) = cell<triangle::N()>({(j+1)*n + i+1, (j+1)*n + i, j*n + i + 1}) ; 
		}
	}
	return m;
}
