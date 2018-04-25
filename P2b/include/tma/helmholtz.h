#ifndef TMA_HELMHOLTZ_H
#define TMA_HELMHOLTZ_H

namespace tma
{

void Assemble(const distributedMesh<triangle>& dm, real (*f)(const point<2>&), uint rank, DistributedSparseMatrix& MassMatrix, DistributedSparseMatrix& StiffnessMatrix, DistributedVector& ForceVector)
{
	triangle T;
	for(uint c = 0; c < dm.ncells(); c++) {
		if(dm.cellOwner(c) != rank) continue;
		T = dm.physicalCell(c);
		
	}
	
}

} // namespace tma

#endif // TMA_HELMHOLTZ_H
