#ifndef TMA_UHEAT_H
#define TMA_UHEAT_H

namespace tma
{

point<2> perpendicular(const point<2>& p)
{
	return point<2>({-p(1), p(0)});
}

void AssembleUnsteadyHeatProblem(const distributedMesh<triangle>& dm, uint rank, DistributedSparseMatrix& StiffnessMatrix, DistributedSparseMatrix& MassMatrix, DistributedSparseMatrix& MixMatrix)
{
	triangle T; real Tsize, mm;
	uint globJ;
	basefunctions<triangle> bf;
	TWBQuadrature quad;
	
	#ifdef OPENMP
	#pragma omp parallel for private(T, Tsize, globJ, bf, quad) schedule(guided)
	#endif
	for(uint v = 0; v < dm.nverts(); v++) {		// go through all vertices
		if(dm.vertOwner(v) != rank) continue;
		//if(dm.geom()(v)(1) == 0. || dm.geom()(v)(1) == 1.) {	// boundary conditions
		//	StiffnessMatrix.set(v, v, 1.);
		//	continue;
		//}
		for(uint c = 0; c < dm.ncells(); c++) {	// find cells the vertex belongs to
			if(dm.topo()(c).contains(v)) {
				T = dm.physicalCell(c);
				Tsize = T.area();
				uint i = dm.topo()(c).find(v);	// local cell index i
				for(uint j = 0; j < 3; j++) {	// find other vertices in cell
					globJ = dm.topo()(c)(j);
					StiffnessMatrix.set(v, globJ, StiffnessMatrix(v, globJ) + (T((i+1)%3)-T((i+2)%3))*(T((j+1)%3)-T((j+2)%3))/4./Tsize);
					MassMatrix.set(v, globJ, MassMatrix(v, globJ) + quad.Integral(T, [&bf, i, j, &T] (const point<2>& x) { return bf.phi(i, T, x)*bf.phi(j, T, x); }));
					
					if((dm.geom()(v)(0) == 0. || dm.geom()(v)(0) == 1. || dm.geom()(v)(1) == 0. || dm.geom()(v)(1) == 1.)	// if cell pair lies on the edge
						&& (dm.geom()(globJ)(0) == 0. || dm.geom()(globJ)(0) == 1. || dm.geom()(globJ)(1) == 0. || dm.geom()(globJ)(1) == 1.))
					{
						mm = 0;
						for(uint k = 0; k < 3; k++) 
							mm += perpendicular(T((j+1)%3) - T((j+2)%3))*perpendicular(T((k+1)%3) - T((k+2)%3)) 
											* bf.phi(i, T, (1./2.)*(T((k+1)%3) + T((k+2)%3))) 	*	 (T((k+1)%3) - T((k+2)%3)).norm();
						MixMatrix.set(v, globJ, MixMatrix(v, globJ) + mm);
					}
				}
			}
		}		
	}	
}

void CalcU_tilde(const distributedMesh<triangle>& dm, const std::function<real(const point<2>&, real)>& u_tildef, const real& t, uint rank, DistributedVector& u_tilde)
{
	triangle T;
	basefunctions<triangle> bf; TWBQuadrature quad;
	#ifdef OPENMP
	#pragma omp parallel for private(T, bf, quad) schedule(guided)
	#endif
	for(uint v = 0; v < dm.nverts(); v++) {
		if(dm.vertOwner(v) != rank) continue;
		u_tilde[v] = 0;
		for(uint c = 0; c < dm.ncells(); c++) {	// find cells the vertex belongs to
			if(dm.topo()(c).contains(v)) {
				T = dm.physicalCell(c);
				uint i = dm.topo()(c).find(v);
				u_tilde[v] += quad.Integral(T, [&bf, &u_tildef, i, &T, &t] (const point<2>& x) { return u_tildef(x, t)*bf.phi(i, T, x);  });
			}
		}	
	}
}

DistributedVector BackwardsEuler(const distributedMesh<triangle>& dm, const DistributedSparseMatrix& M, const DistributedSparseMatrix& K, const DistributedSparseMatrix& S,
					const std::function<real(const point<2>&, real)>& u_tildef, const real& nu, const real& h, const real& t_0, uint rank, uint nit)
{
	DistributedVector u(dm, rank), u_tilde(dm, rank), Mu(dm, rank), Su_tilde(dm, rank);
	DistributedSparseMatrix MhnuK(M + h*nu*K);
	real t = t_0;
	CalcU_tilde(dm, u_tildef, t, rank, u_tilde);
	Sync(dm, u_tilde, rank);
	u = u_tilde;
	
	std::cout << rank << ": u_0 = " << u << std::endl;
	
	for(uint i = 0; i < nit; i++)
	{
		t+= h;
		CalcU_tilde(dm, u_tildef, t, rank, u_tilde); Sync(dm, u_tilde, rank);
		
		DistributedMatrixVectorProduct(dm, M, u, Mu); Sync(dm, Mu, rank);
		DistributedMatrixVectorProduct(dm, S, u_tilde, Su_tilde); Sync(dm, Mu, rank);
		u = ConjugateGradient(dm, MhnuK, Mu + Su_tilde, rank, 1e-4, 500);
	}
	
	return u;
}



} // namespace tma


#endif // TMA_UHEAT_H
