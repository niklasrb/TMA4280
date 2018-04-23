#ifndef TMA_MATRIX_H
#define TMA_MATRIX_H

#include <cmath>
#include <algorithm>

namespace tma
{
	
class SparseMatrix
{
protected:
	uint m, n, nnz;
	std::vector<uint> IA, JA;
	std::vector<real> A;
	
public:
	SparseMatrix(uint m, uint n, double** matrix) : m(m), n(n)
	{
		nnz = 0;
		IA.push_back(nnz);
		for(uint i = 0; i < m; i++) {
			for(uint j = 0; j < n; j++) {
				if(matrix[i][j] != 0) {
					A.push_back(matrix[i][j]);
					JA.push_back(j);
					nnz++;
				}
			}
			IA.push_back(nnz);
		} 
	}
	SparseMatrix(uint m, uint n) : m(m), n(n)
	{
		nnz = 0;
		IA.resize(m+1, 0);
	}
	
	SparseMatrix() : m(0), n(0)
	{
		
	}
	
	friend std::ostream& operator <<(std::ostream& os,const SparseMatrix& sm)
	{
		//os << "A: ("; for(unsigned int i = 0; i < sm.A.size(); i++) os << sm.A[i] << (i < sm.A.size()-1 ? ", " : "");  os << ")" << std::endl;
		//os << "IA: ("; for(unsigned int i = 0; i < sm.IA.size(); i++) os << sm.IA[i] << (i < sm.IA.size()-1 ? ", " : ""); os << ")" << std::endl;
		//os << "JA: ("; for(unsigned int i = 0; i < sm.JA.size(); i++) os << sm.JA[i] << (i < sm.JA.size()-1 ? ", " : ""); os << ")" << std::endl;
		os << "[ ";
		for(uint i = 0; i < sm.m; i++) {
			os << "[ ";
			for(uint j = 0; j < sm.n; j++)
				os << sm(i, j) << (j < sm.n - 1 ? ", " : "]");
			if(i < sm.m-1) os << "," << std::endl;
		}
		return os << "]";
	}
	
	void set(uint i, uint j, const double& v)
	{
		assert(0 <= i && i < m && 0 <= j && j < n);
		if((*this)(i, j) != 0) {	// element exists, find it
			int c = IA.at(i);
			while(JA.at(c) != j) c++;
			if(v == 0) {	// delete it
				A.erase(A.begin()+c);
				JA.erase(JA.begin()+c);
				nnz--;
				for(uint k = i+1; k < m+1; k++) IA.at(k)--;
			}
			else {			//change it
				A[c] = v;
			}
		}
		else {		// element does not exist
			if(v == 0) return; // no need to do anything
			// create it
			uint c = IA.at(i);
			while(c < JA.size() && c < IA.at(i+1) && JA.at(c) < j) c++;
			A.insert(A.begin() + c , v);
			JA.insert(JA.begin()+ c , j);
			nnz++;
			for(uint k = i+1; k < m+1; k++) IA.at(k)++;
		}
	}
	
	real operator ()(uint i, uint j) const
	{
		assert(0 <= i && i < m && 0 <= j && j < n);
		for(uint c = IA.at(i); c < IA.at(i+1); c++)
			if(JA.at(c) == j) return A[c];
		return 0;
	}
	
	friend SparseMatrix operator +(const SparseMatrix& sm1, const SparseMatrix& sm2)
	{
		assert(sm1.m == sm2.m);
		assert(sm1.n == sm2.n);
		SparseMatrix res(sm1.m, sm1.n);
		real d;
		for(uint i = 0; i < sm1.m; i++) {
			for(uint j = 0; j < sm1.n; j++) {
				d = sm1(i, j) + sm2(i, j);
				if(d != 0) res.set(i, j, d);
			}				
		}
		return res;
	}
	
	friend SparseMatrix operator -(const SparseMatrix& sm) 
	{
		SparseMatrix res(sm);
		for(uint i = 0; i < res.A.size(); i++)
			res.A[i] = -res.A[i];
		return res;
	}
	
	friend SparseMatrix operator *(const SparseMatrix& sm1, const SparseMatrix& sm2)
	{
		assert(sm1.n == sm2.m);
		SparseMatrix res(sm1.m, sm2.n);
		real d;
		for(uint i = 0; i < sm1.m; i++) {
			for(uint j = 0; j < sm2.n; j++) {
				d = 0;
				for(uint k = 0; k < sm1.n; k++)
					d += sm1(i, k)*sm2(k, j);
				if(d!=0)
					res.set(i, j, d);
			}
		}
	}
};

class DistributedSparseMatrix
{
protected:
	SparseMatrix diag_, offdiag_;
	map<uint, uint> globalToLocalRow_;
	std::pair<uint, uint> range_;
	uint m, n;

public:
	template<class T, uint D>
	DistributedSparseMatrix(const distributedMesh<T, D>& dM, uint n, uint rank) : n(n)
	{
		uint r = 0;
		for(uint v = 0; v < dM.nverts(); v++) {
			if(dM.vertOwner(v) == rank) 
				globalToLocalRow_[v] = r++;
		}
		diag_ = SparseMatrix(r, n);
		offdiag_ = SparseMatrix(r, n);
		m = r;
		bool range = true;
		for(uint i = 0; i < r-1; i++) {
			if(LocalToGlobalRowIndex(i)+1 != LocalToGlobalRowIndex(i+1))
				range = false;
		}
		if(range)
			range_ = std::make_pair(LocalToGlobalRowIndex(0), LocalToGlobalRowIndex(r-1)+1);
		else
			range_ = std::make_pair(std::nan(""), std::nan(""));
	}
	
protected:
	uint LocalToGlobalRowIndex(uint i) const 
	{
		for(auto it = globalToLocalRow_.begin(); it != globalToLocalRow_.end(); ++it)
			if(it->second == i) return it->first;
		std::cout << "couldn't find " << i << " here :(" << std::endl;
		assert(false);
	}

	// Takes local indices and checks wether they are diagonal in the complete matrix
	bool IsInDiagonalBlock(uint i, uint j) const 
	{
		assert(range_.first == range_.first);
		if(range_.first != range_.first) 
			return true;
		return (range_.first <= i && range_.first <= j && i < range_.second && j < range_.second);
	}

public:	
	real operator()(uint i, uint j) const
	{
		assert(globalToLocalRow_.find(i) != globalToLocalRow_.end() && j < n);
		if(IsInDiagonalBlock(i, j))
			return diag_(globalToLocalRow_.at(i), j);
		else
			return offdiag_(globalToLocalRow_.at(i), j);
	}
	
	void set(uint i, uint j, real v)
	{
		assert(globalToLocalRow_.find(i) != globalToLocalRow_.end() && j < n);
		if(IsInDiagonalBlock(i, j))
			diag_.set(globalToLocalRow_.at(i), j, v);
		else
			offdiag_.set(globalToLocalRow_.at(i), j, v);
	}
	
	void set(std::pair<uint, uint> index, real v) { set(index.first, index.second, v); }
	real operator()(std::pair<uint, uint> index) const { return (*this)(index.first, index.second); }
	
	friend std::ostream& operator <<(std::ostream& os, const DistributedSparseMatrix& dsm) 
	{
		return os << dsm.diag_ + dsm.offdiag_;
	}
	
	std::pair<uint, uint> RowRange() const { return range_; }
	uint col() const { return n; }

};
	
	
} /* namspace tma */



#endif /* TMA_MATRIX_H */
