#ifndef TMA_MESH_H_
#define TMA_MESH_H_

#include <tma/types.h>

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <vector>

namespace tma
{

/*
 * This is a simple topology class to store cell-vertices connectivities in
 * an array and return the number of entities.
 * To make it easier, a template argument is used: it defines the cell type.
 */
template <class T>
struct topology
{
  // The constructor requires the number of cells and vertices
  topology(uint ncells, uint nverts) :
    ncells_(ncells),
    nverts_(nverts),
    cv_(ncells ? new uint[ncells * T::num(0)]() : NULL), // Need to make sure it is non-zero
    offset_(0)
  {
  }

  topology(const topology& t) :
	ncells_(t.ncells_),
	nverts_(t.nverts_),
	cv_(t.ncells_ ? new uint[t.ncells_ * T::num(0)]() : NULL), // Need to make sure it is non-zero
    offset_(0)
    {
		for(uint i = 0; i < ncells_*T::num(0); i++)
			cv_[i] = t.cv_[i];
	}

  ~topology()
  {
    delete [] cv_;
  }

  // Return the topological dimension
  uint dim() const { return T::dim(); }

  // Return the number of cells in the mesh
  uint ncells() const { return ncells_; }

  // Return the number of vertices in the mesh
  uint nverts() const { return nverts_; }
  
  // Return the number of vertices in a cell
  uint nnum() const { return T::num(0); }

  // Accessor for the vertices of the i-th cell
  uint * operator()(uint i)
  {
    return cv_ + i * T::num(0);
  }
 

  // Display the connectivities
  void dump() const
  {
    uint const * cv = cv_;
    for (uint c = 0; c < ncells_; ++c)
    {
      std::cout << std::setw(4) << c << ":";
      for (uint v = 0; v < T::num(0); ++v, ++cv)
      {
        std::cout << std::setw(4) << *cv;
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

private:

  uint const ncells_;
  uint const nverts_;
  uint * const cv_;
  uint offset_; // Use in parallel

};



/*
 * This is a simple geometry class to store vertex coordinates in an array.
 */
template <uint D>
struct geometry
{
  geometry(uint nverts) :
    nverts_(nverts),
    vx_(nverts ? new real[nverts*D]() : NULL)
  {

  }
  
  geometry(const geometry& gm) : 
	nverts_(gm.nverts_), vx_(gm.nverts_ ? new real[gm.nverts_*D]() : NULL)
	{
		for(uint i = 0; i < nverts_; i++)
			vx_[i] = gm.vx_[i];
	}

  ~geometry()
  {
    delete [] vx_;
  }

  // Accessor for the coordinates of the i-th vertex
  real * operator()(uint i)
  {
    return vx_ + i * D;
  }

  // Display the coordinates
  void dump() const
  {
    real const * vx = vx_;
    for (uint v = 0; v < nverts_; ++v)
    {
      std::cout << std::setw(4) << v << ":";
      for (uint x = 0; x < D; ++x, ++vx)
      {
        std::cout << std::setw(8) << *vx;
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

private:

  uint const nverts_;
  real * const vx_;

};


template <class T, uint D>
class mesh
{
public:
  /*
   * Now let us define a mesh class combining both ingredients to build a mesh
   * for a given cell type in a given Euclidean space.
   * We use template arguments to define the topology and the geometry and we
   * only need to specify the number of cells and vertices.
   */

  mesh(uint ncells, uint nverts) :
    T_(ncells, nverts),
    G_(nverts)
  {
  }

  topology<T>& topo() { return T_; }
  geometry<D>& geom() { return G_; }

protected:

  topology<T> T_;
  geometry<D> G_;
};

template <class T, uint D>
class distributedMesh : 
	public mesh<T, D>
{
public:
	distributedMesh(uint ncells, uint nverts, uint nprocs) : mesh<T, D>(ncells, nverts), nprocs_(nprocs)
	{
		int nCellsPerProcess = ceil(double(ncells)/nprocs);
		uint r = 0, c = 0;
		while(r < nprocs && c < ncells) {
			cellOwner_.insert(std::make_pair(c , r));
			if(++c % nCellsPerProcess == 0) r++;
		}
	}
	
	distributedMesh(const mesh<T, D>& m, uint nprocs) : mesh<T, D>(m), nprocs_(nprocs)
	{
		int nCellsPerProcess = ceil(double(this->T_.ncells())/nprocs);
		uint r = 0, c = 0;
		while(r < nprocs && c < this->T_.ncells()) {
			cellOwner_.insert(std::make_pair(c , r));
			if(++c % nCellsPerProcess == 0) r++;
		}
		AssignVerticesOwner();
	}
	
	void AssignVerticesOwner()
	{
		uint owner;
		for(uint c = 0; c < this->T_.ncells(); c++) {
			owner = cellOwner_.at(c);
			for(uint i = 0; i < T::num(0); i++) {
				auto v = vertOwner_.find(this->T_(c)[i]);
				if(v == vertOwner_.end()) 
					vertOwner_.insert(std::make_pair(this->T_(c)[i], owner));
				else if(owner < v->second )
					vertOwner_[v->first] = owner;
			}
		}
	}
	
	void dump()
	{	// T_ dump
		for (uint c = 0; c < ncells(); ++c)  
		{
			std::cout << std::setw(4) << c << " owned by " << cellOwner(c) << " :";
			for (int v = 0; v < this->T_.nnum(); ++v)
				std::cout << std::setw(4) << this->T_(c)[v] << (v == this->T_.nnum()-1 ? '\n' : ' ');
		} std::cout << std::endl;
		// G_ dump
		for (uint v = 0; v < nverts(); ++v)
		{
			std::cout << std::setw(4) << v << " owned by " << vertOwner(v) << " :";
			for (uint x = 0; x < D; ++x)
				std::cout << std::setw(8) << this->G_(v)[x] << (x == D-1 ? '\n' : ' ');
		} std::cout << std::endl;
	}
	
	
	uint ncells() const { return this->T_.ncells(); }
	uint nverts() const { return this->T_.nverts(); }
	uint cellOwner(uint cell) const { return cellOwner_.at(cell); }
	uint vertOwner(uint vertex) const { return vertOwner_.at(vertex); }
	
protected:
	std::map<uint, uint> cellOwner_, vertOwner_;
	uint nprocs_;
};

template<typename key, typename value> 
class map : public std::map<key, value>
{
public:
	bool contains(value v) {
		for(auto it = this->begin(); it != this->end(); ++it) {
			if(it->second == v) 
				return true;
		}
		return false;
	}
};

template<typename Z>
class distributedVector
{
protected:
	std::vector<Z> local_, ghosts_;
	map<uint, uint> localToGlobal_;
	std::pair<uint, uint> range_;
	
public:
	template<class T, uint D>
	distributedVector(distributedMesh<T, D>& dM, uint rank)
	{
		for(uint v = 0; v < dM.nverts(); v++) {
			if(dM.vertOwner(v) == rank) {
				local_.push_back(0);
				localToGlobal_[local_.size()-1] = v;
			}
		}
		for(uint c = 0; c < dM.ncells(); c++) {
			if(dM.cellOwner(c) == rank) {
				for(uint i = 0; i < T::num(0); i++) {
					uint v = dM.topo()(c)[i];
					if(!localToGlobal_.contains(v)) {
						ghosts_.push_back(0);
						localToGlobal_[local_.size() + ghosts_.size() -1] = v;
					}
				}
			}
		}
	}
	
	Z& operator[](uint i) {
		assert(i < local_.size() + ghosts_.size());
		if(i < local_.size()) return local_[i];
		return ghosts_[i-local_.size()];
	}
	
	
	uint size() const { return local_.size(); }
	uint sizeTotal() const { return local_.size() + ghosts_.size(); }
	uint LocalToGlobalIndex(uint i) const { return localToGlobal_.at(i); }
	uint GlobalToLocalIndex(uint i) const 
	{
		for(uint j = 0; j < local_.size() + ghosts_.size(); j++)
			if(localToGlobal_.at(j) == i) return j;
		assert(false);
	}
	std::pair<uint, uint> Range() const { return range_; }
};

} /* namespace tma */

#endif /* TMA_MESH_H_ */
