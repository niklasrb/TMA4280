#ifndef TMA_HELMHOLTZ_H
#define TMA_HELMHOLTZ_H

namespace tma
{

void AssembleHelmholtz(const distributedMesh<triangle>& dm, real (*f)(const point<2>&), uint rank, DistributedSparseMatrix& StiffnessMatrix, DistributedSparseMatrix& MassMatrix, DistributedVector& ForceVector)
{
	triangle T; real Tsize;
	uint globJ;
	basefunctions<triangle> bf;
	TWBQuadrature quad;
	
	#ifdef OPENMP
	#pragma omp parallel for private(T, Tsize, globJ, bf, quad) schedule(guided)
	#endif
	for(uint v = 0; v < dm.nverts(); v++) {		// go through all vertices
		if(dm.vertOwner(v) != rank) continue;
		if(dm.geom()(v)(1) == 0. || dm.geom()(v)(1) == 1.) {	// boundary conditions
			StiffnessMatrix.set(v, v, 1.);
			continue;
		}
		for(uint c = 0; c < dm.ncells(); c++) {	// find cells the vertex belongs to
			if(dm.topo()(c).contains(v)) {
				T = dm.physicalCell(c);
				Tsize = T.area();
				uint i = dm.topo()(c).find(v);
				for(uint j = 0; j < 3; j++) {	// find other vertices in cell
					globJ = dm.topo()(c)(j);
					StiffnessMatrix.set(v, globJ, StiffnessMatrix(v, globJ) + (T((i+1)%3)-T((i+2)%3))*(T((j+1)%3)-T((j+2)%3))/4./Tsize);
					MassMatrix.set(v, globJ, MassMatrix(v, globJ) + quad.Integral(T, [&bf, i, j, &T] (const point<2>& x) { return bf.phi(i, T, x)*bf.phi(j, T, x); }));
				}
				ForceVector[v] += quad.Integral(T, [&bf, &f, i, &T] (const point<2>& x) { return f(x)*bf.phi(i, T, x);  });
			}
		}	
		
	}
	
}

} // namespace tma

#endif // TMA_HELMHOLTZ_H
