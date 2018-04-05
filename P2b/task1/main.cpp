
#include "../include/tma.h"

using namespace tma;

int main(int argc, char**argv)
{
	// Test mesh constructor
	int N = 15, P =4;
	
	mesh<interval, 1> m(N, N+1);
	for(int i = 0; i < N; i++ ) {
		m.topo()(i)[0] = i;
		m.topo()(i)[1] = i+1;
	}
	for(int i = 0; i < N+1; i++)
		m.geom()(i)[0] = double(i)/N;
	
	m.topo().dump();
	m.geom().dump();
	distributedMesh<interval, 1> dm(m, P);
	dm.dump();
	
	distributedVector<real> dV0(dm, 0);
	std::cout << "distributed Vector rank 0" << std::endl;
	for(uint i = 0; i < dV0.size(); i++)
		std::cout << "\tlocal: " << i << " <-> global " << dV0.LocalToGlobalIndex(i) << std::endl;
	std::cout << "\tghosts: " << std::endl;
	for(uint i = dV0.size(); i < dV0.sizeTotal(); i++)
		std::cout << "\tlocal: " << i << " <-> global " << dV0.LocalToGlobalIndex(i) << std::endl;
	
	distributedVector<real> dV1(dm, 3);
	std::cout << "distributed Vector rank 3" << std::endl;
	for(uint i = 0; i < dV1.size(); i++)
		std::cout << "\tlocal: " << i << " <-> global " << dV1.LocalToGlobalIndex(i) << std::endl;
	std::cout << "\tghosts: " << std::endl;
	for(uint i = dV1.size(); i < dV1.sizeTotal(); i++)
		std::cout << "\tlocal: " << i << " <-> global " << dV1.LocalToGlobalIndex(i) << std::endl;
	
	return 0;
}
