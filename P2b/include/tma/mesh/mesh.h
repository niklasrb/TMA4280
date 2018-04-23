#ifndef TMA_MESH_H_
#define TMA_MESH_H_

#include <tma/types.h>

#include <iostream>
#include <iomanip>
#include <cmath>


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
		for(uint i = 0; i < D*nverts_; i++)
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
protected:
	map<uint, uint> cellOwner_, vertOwner_;
	map<uint, vector<uint> > ghosts_;
	uint nprocs_;
public:
	/*
	 *  Create empty mesh and distribute cells among processes
	 *  Since cells don't contain vertices yet, they have to be assigned to owner later
	 *
	distributedMesh(uint ncells, uint nverts, uint nprocs) : mesh<T, D>(ncells, nverts), nprocs_(nprocs)	
	{	
		int nCellsPerProcess = ceil(double(ncells)/nprocs);
		uint r = 0, c = 0;
		while(r < nprocs && c < ncells) {
			cellOwner_.insert(std::make_pair(c , r));
			if(++c % nCellsPerProcess == 0) r++;
		}
	}*/
	/*
	 * Use an exisiting mesh to distribute cells and vertices among processes
	 */
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
	
	/*
	 * Assigns vertices to their owner according to cell ownership
	 */
	void AssignVerticesOwner()
	{
		//std::cout << "AssignVerticesOwner()" << std::endl;
		vertOwner_.clear();
		uint owner;
		for(uint c = 0; c < this->T_.ncells(); c++) {
			owner = cellOwner_.at(c);					// owner of this cell
			//std::cout << "cell " << c << " with owner " << owner << std::endl;
			for(uint i = 0; i < T::num(0); i++) {
				auto vertex = this->T_(c)[i];
				auto vit = vertOwner_.find(vertex);
				if(vit == vertOwner_.end()) {				// vertex has no owner so far
					//std::cout << "vertex " << vertex << " has no owner so far, so it will be " << owner << std::endl;
					vertOwner_.insert(std::make_pair(this->T_(c)[i], owner));
				}
				else if (vit->second != owner) {									// vertex has other owner
					auto git = ghosts_.find(vertex);
					uint ghost = owner;
					if(owner < vit->second ) {				// this owner has lower rank 
						vertOwner_[vertex] = owner;
						ghost = vit->second;	
					}
					if(git == ghosts_.end())				// add the secondary owner to list
						ghosts_.insert(std::make_pair(vertex, vector<uint>()));
					ghosts_[vertex].push_back(ghost);	
					//std::cout << "vertex " << vertex << " has other owner " << vit->second  << " so far, so " << ghost  << " will be a ghost " << std::endl; 
				}
			}
		}
		// If you own a vertex in a cell, you should be able to see the whole cell
		for(uint c = 0; c < this->T_.ncells(); c++) {
			for(uint i = 0; i < T::num(0); i++) {
				auto vertex = this->T_(c)[i];		
				for(uint j = 0; j < T::num(0); j++) {
					if(i == j)
						continue;
					if(vertOwner_[vertex] == vertOwner_[this->T_(c)[j]])
						continue;
					if(!ghosts_[this->T_(c)[j]].contains(vertOwner_[vertex]))
						ghosts_[this->T_(c)[j]].push_back(vertOwner_[vertex]);
				}
			}
		}
	}
	
	void dump()
	{	// T_ dump
		for (uint c = 0; c < ncells(); ++c)  
		{
			std::cout << std::setw(4) << c << " owned by " << cellOwner(c) << " :";
			for (uint v = 0; v < this->T_.nnum(); ++v)
				std::cout << std::setw(4) << this->T_(c)[v] << (v == this->T_.nnum()-1 ? '\n' : ' ');
		} std::cout << std::endl;
		// G_ dump
		for (uint v = 0; v < nverts(); ++v)
		{
			std::cout << std::setw(4) << v << " owned by " << vertOwner(v);
			 if(ghosts_.find(v) != ghosts_.end())
				std::cout <<"(shared with " << ghosts_[v] << ")";
			std::cout << ": ";
			for (uint x = 0; x < D; ++x)
				std::cout << std::setw(8) << this->G_(v)[x] << (x == D-1 ? '\n' : ' ');
		} std::cout << std::endl;
	}
	
	
	uint ncells() const { return this->T_.ncells(); }
	uint nverts() const { return this->T_.nverts(); }
	uint cellOwner(uint cell) const { return cellOwner_.at(cell); }
	uint vertOwner(uint vertex) const { return vertOwner_.at(vertex); }
	bool IsGhost(uint vertex) const { return ghosts_.find(vertex) != ghosts_.end(); }
	vector<uint> secondaryOwners(uint vertex) const { auto vit = ghosts_.find(vertex); return (vit == ghosts_.end() ? vector<uint>() : vit->second ); }

	
};


} /* namespace tma */

#endif /* TMA_MESH_H_ */
