#ifndef TMA_VECTOR_H
#define TMA_VECTOR_H

namespace tma
{

	
class DistributedVector
{
protected:
	std::vector<real> local_, ghosts_;
	map<uint, int> globalToLocal_;
	std::pair<uint, uint> range_;
	
	uint LocalToGlobalIndex(int i) const 
	{
		for(auto it = globalToLocal_.begin(); it != globalToLocal_.end(); ++it)
			if(it->second == i) return it->first;
		assert(false);
	}
	
public:
	template<class T, uint D>
	DistributedVector(const distributedMesh<T, D>& dM, uint rank)	// use mesh to initiate vector
	{
		for(uint v = 0; v < dM.nverts(); v++) {		// find owned vertices
			if(dM.vertOwner(v) == rank) {
				local_.push_back(0);
				globalToLocal_[v] = local_.size()-1;
			} else if(dM.IsGhost(v)) {				// and ghosts for this instance
				if(dM.secondaryOwners(v).contains(rank)) {
					ghosts_.push_back(v);
					globalToLocal_[v] = -(ghosts_.size()-1)-1;
				}
			}
		}
		bool range = true;
		for(uint i = 0; i < local_.size()-1; i++) {
			if(LocalToGlobalIndex(i)+1 != LocalToGlobalIndex(i+1))
				range = false;
		}
		if(range)
			range_ = std::make_pair(LocalToGlobalIndex(0) , LocalToGlobalIndex(local_.size()-1)+1);
		else
			range_ = std::make_pair(std::nan(""), std::nan(""));
		
	}
	
	// access elements
	real& operator[](uint i) {		
		assert(globalToLocal_.find(i) != globalToLocal_.end());
		assert(globalToLocal_[i] >= 0);
		return local_.at(globalToLocal_[i]);
	}
	
	real operator[](uint i) const {		
		assert(globalToLocal_.find(i) != globalToLocal_.end());
		if(globalToLocal_.at(i) >= 0) return local_.at(globalToLocal_.at(i));
		return ghosts_.at(-globalToLocal_.at(i)-1);
	}
	
	real at(uint i) const {
		return (*this)[i];
	}
	
	// Lets you update ghost entries
	void update(uint i, real v) 
	{ 
		assert(globalToLocal_.find(i)!= globalToLocal_.end() && globalToLocal_[i] < 0);
		ghosts_.at(-globalToLocal_[i]-1) = v;
	}
	
	// returns number of local elements
	uint size() const { return local_.size(); }
	
	// returns number of total elements, local + ghosts
	uint sizeTotal() const { return local_.size() + ghosts_.size(); }
	
	// shows the global range of this distributed vector, if exists
	std::pair<uint, uint> Range() const { return range_; }
	
	friend std::ostream& operator <<(std::ostream& os, const DistributedVector& v)
	{
		os << "[";
		for(uint i = 0; i < v.local_.size(); i++)
			os << v.LocalToGlobalIndex(i) << ":" << v[v.LocalToGlobalIndex(i)] << (i < v.local_.size()-1 ? "\t" : (v.ghosts_.size() > 0 ? " ; " : "]"));
		for(int i = 1; i <= v.ghosts_.size(); i++)
			os << v.LocalToGlobalIndex(-i) << ":" << v[v.LocalToGlobalIndex(-i)] << (i < v.ghosts_.size()-1 ? "\t" : "]");
		return os;
	}
	
	friend DistributedVector operator +(const DistributedVector& v, const DistributedVector& w)
	{
		assert(v.local_.size() == w.local_.size());
		assert(v.ghosts_.size() == w.ghosts_.size());
		assert(v.globalToLocal_ == w.globalToLocal_);
		DistributedVector u(v);
		for(uint i = 0; i < w.local_.size(); i++)
			u.local_[i] += w.local_[i];
		for(uint i = 0; i < w.ghosts_.size(); i++)
			u.ghosts_[i] += w.ghosts_[i];
		return u;
	}
	
	friend DistributedVector operator -(const DistributedVector& v)
	{
		DistributedVector u(v);
		for(uint i = 0; i < u.local_.size(); i++)
			u.local_[i] *= -1;
		for(uint i = 0; i < u.ghosts_.size(); i++)
			u.ghosts_[i] *= -1;
		return u;
	}
	
	real norm() const 
	{ 
		real s = 0; 
		for(uint i = 0; i < local_.size(); i++) 
			s += pow(local_[i],2); 
		return sqrt(s);  
	}
	
	void dump() const
	{
		std::cout << (*this) << "-(" << range_.first << ", " << range_.second << ")" << std::endl;
	}
};
	
} /* namespace tma */

#endif /* TMA_VECTOR_H */
