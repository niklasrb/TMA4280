#ifndef TMA_MESH_H_
#define TMA_MESH_H_

#include <tma/types.h>

#include <iostream>
#include <iomanip>
#include <cmath>


namespace tma
{
template <uint D>
struct geometry
{
protected:
	vector<point<D> > vertices_;
public:
  geometry(uint nverts) 
    {
		vertices_.resize(nverts);
	}

  point<D>& operator()(uint i) {  return vertices_.at(i);  }
  point<D> operator()(uint i) const {  return vertices_.at(i);  }
  
  uint nverts() const {  return vertices_.size();  }

  // Display the coordinates
  void dump() const
  {
     for(uint i = 0; i < vertices_.size(); i++)
		std::cout << i << ": " << vertices_.at(i) << std::endl;
  }

};

template <class T>
struct topology
{
protected:
	vector<cell<T::N()> > cells_;
public:
  // The constructor requires the number of cells and vertices
	topology(uint ncells) {	cells_.resize(ncells); }

  // Return the number of cells in the mesh
  uint ncells() const { return cells_.size(); }

  // Accessor for the vertices of the i-th cell
  cell<T::N()>& operator()(uint i) { return cells_.at(i); }
  cell<T::N()> operator()(uint i) const { return cells_.at(i); }

  // Display the connectivities
  void dump() const
  {
		for(auto it = cells_.begin(); it != cells_.end(); ++it)
			std::cout << (*it) << std::endl;
  }

};


template<class T>
class mesh
{
protected:
	geometry<T::D()> geom_;
	topology<T> topo_;
	
public:
	uint nverts() const { return geom_.nverts(); }
	uint ncells() const { return topo_.ncells(); }
	
	mesh(uint ncells, uint nverts) : geom_(nverts), topo_(ncells) { std::cout << ncells << " - " << nverts << std::endl;}
	
	geometry<T::D()>& geom() { return geom_; }
	geometry<T::D()> geom() const { return geom_; }
	topology<T>& topo() { return topo_; }
	topology<T> topo() const { return topo_; }
	
	
	T physicalCell(uint c) const { return physicalCell(topo_(c)); }
	
	T physicalCell(const cell<T::N()>& c) const
	{
		vector<point<T::D()> > p;
		for(uint i = 0; i < c.size(); i++) p.push_back( geom_(c(i)));
		return T(p);
	}
	
	void dump() const
	{
		std::cout << "geometry: " << std::endl;
		geom_.dump();
		std::cout << "topology: " << std::endl;
		topo_.dump();
	}
};


template <class T>
class distributedMesh : 
	public mesh<T>
{
protected:
	map<uint, uint> cellOwner_, vertOwner_;
	map<uint, vector<uint> > ghosts_;
	uint nprocs_;
public:
	/*
	 * Use an exisiting mesh to distribute cells and vertices among processes
	 */
	distributedMesh(const mesh<T>& m, uint nprocs) : mesh<T>(m), nprocs_(nprocs)
	{
		int nCellsPerProcess = ceil(double(this->ncells())/nprocs);
		uint r = 0, c = 0;
		while(r < nprocs && c < this->ncells()) {
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
		for(uint c = 0; c < this->ncells(); c++) {
			owner = cellOwner_.at(c);					// owner of this cell
			//std::cout << "cell " << c << " with owner " << owner << std::endl;
			for(uint i = 0; i < this->topo()(c).size(); i++) {
				auto vertex = this->topo()(c)(i);
				auto vit = vertOwner_.find(vertex);
				if(vit == vertOwner_.end()) {				// vertex has no owner so far
					//std::cout << "vertex " << vertex << " has no owner so far, so it will be " << owner << std::endl;
					vertOwner_.insert(std::make_pair(this->topo()(c)(i), owner));
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
					if(!ghosts_[vertex].contains(ghost))
						ghosts_[vertex].push_back(ghost);	
					//std::cout << "vertex " << vertex << " has other owner " << vit->second  << " so far, so " << ghost  << " will be a ghost " << std::endl; 
				}
			}
		}
		// If you own a vertex in a cell, you should be able to see the whole cell
		for(uint c = 0; c < this->ncells(); c++) {
			for(uint i = 0; i < this->topo()(c).size(); i++) {
				auto vertex = this->topo()(c)(i);	
				auto owner = vertOwner_[vertex];
				for(uint j = 0; j < this->topo()(c).size(); j++) {
					if(i == j)
						continue;
					if(owner == vertOwner_[this->topo()(c)(j)])
						continue;
					if(!ghosts_[this->topo()(c)(j)].contains(owner))
						ghosts_[this->topo()(c)(j)].push_back(owner);
				}
			}
		}
	}
	
	void dump() const
	{	// T_ dump
		for (uint c = 0; c < this->ncells(); ++c)  
		{
			std::cout << std::setw(4) << c << " owned by " << cellOwner(c) << " :" << this->topo()(c) << std::endl;
		} 
		// G_ dump
		for (uint v = 0; v < this->nverts(); ++v)
		{
			std::cout << std::setw(4) << v << " owned by " << vertOwner(v);
			 if(ghosts_.find(v) != ghosts_.end())
				std::cout <<"(shared with " << ghosts_.at(v) << ")";
			std::cout << ": " << this->geom()(v) << std::endl;
		}
	}
	
	uint cellOwner(uint cell) const { return cellOwner_.at(cell); }
	uint vertOwner(uint vertex) const { return vertOwner_.at(vertex); }
	bool IsGhost(uint vertex) const { return ghosts_.find(vertex) != ghosts_.end(); }
	vector<uint> secondaryOwners(uint vertex) const { auto vit = ghosts_.find(vertex); return (vit == ghosts_.end() ? vector<uint>() : vit->second ); }

	
};


} /* namespace tma */

#endif /* TMA_MESH_H_ */
