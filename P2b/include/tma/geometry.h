#ifndef TMA_GEOMETRY_H
#define TMA_GEOMETRY_H

#include "types.h"		


namespace tma
{

template<uint D>
struct point
{
protected:
	vector<real> x;
public:
	real& operator ()(uint i) { return x.at(i); }
	real operator()(uint i) const { return x.at(i); }
	point(const vector<real>& p) : x(p) { assert(x.size() == D); }
	point(real a) : x({a}) { assert(x.size() == 1); }
	point()  { x.resize(D); }
	
	friend std::ostream& operator <<(std::ostream& os, const point& I)  
	{  
		if(D > 1)
			os << "[ ";
		for(uint i = 0; i < D; i++)
			os << I(i) << (i < D-1 ? ", " : (D > 1 ? "]" : ""));
		return os;
	}
	
	friend point operator +(const point& p1, const point& p2) { point p(p1); for(uint i = 0; i < D; i++) p(i) += p2(i); return p; } 
	friend point operator -(const point& p1, const point& p2) { point p(p1); for(uint i = 0; i < D; i++) p(i) -= p2(i); return p; } 
	friend point operator *(const real& s, const point& p1) { point p(p1); for(uint i = 0; i < D; i++) p(i) *= s; return p; } 
	friend bool operator ==(const point& p, const point& q) { for(uint i = 0; i < p.x.size(); i++) if(p(i) != q(i)) return false; return true; }
	//point& operator =(std::initializer_list<real> il) {
	friend real operator *(const point& p1, const point& p2) { real s = 0; for(uint i = 0; i < D; i++) s+= p1(i)*p2(i); return s; }
	real norm() const { return std::sqrt((*this)*(*this)); } 
};


struct interval
{
	point<1> a, b;
	constexpr static uint D() { return 1; }
	constexpr static uint N() { return 2; }
	interval(const point<1>& a, const point<1>& b) : a(a), b(b) { }
	interval(const vector<point<1> >& v) { a = v.at(0); b = v.at(1); }
	interval() : a({0}), b({1}) {}
	point<1>& operator ()(uint i) { assert(i < 2); if(i == 0) return a; if(i == 1) return b; }
	point<1> operator ()(uint i) const { assert(i < 2); if(i == 0) return a; if(i == 1) return b; }
	friend std::ostream& operator <<(std::ostream& os, const interval& I) { return os << "(" << I.a << ", " << I.b << ")"; }
	real area() const { return std::abs(b(0) - a(0)); }
};

struct triangle
{
private:
	void sortVertices() 
	{ 
		vector<real> m({(c - b).norm(), (c - a).norm(), (a - b).norm()}); int i = m.maxIndex();  
		vector<point<2> > v({a, b, c}); a = v[i]; b = v[(i+1)%3]; c = v[(i+2)%3]; 
	}
public:
	point<2> a, b, c;
	constexpr static uint D() { return 2; }
	constexpr static uint N() { return 3; }
	triangle(const point<2>& a, const point<2>& b, const point<2>& c) : a(a), b(b), c(c) { sortVertices();}
	triangle(const vector<point<2> >& v) { a = v.at(0); b = v.at(1); c = v.at(2); sortVertices(); }
	triangle() : a({0, 0}), b({1,0}), c({0,1}) {}
	point<2>& operator ()(uint i) { assert(i < 3); if(i ==0) return a; if(i == 1) return b; return c; } 
	point<2> operator ()(uint i) const { assert(i < 3); if(i ==0) return a; if(i == 1) return b; if(i == 2) return c; } 
	friend std::ostream& operator <<(std::ostream& os, const triangle& T) { return os << "( " << T.a << ", " << T.b << ", " << T.c << ")"; }
	real area() const { return (b-a).norm() * (c-a).norm() / 2.; }
};

template<uint N>
struct cell
{
protected:
	vector<uint> vertexIndices_;
public:
	cell(const vector<uint>& vI) : vertexIndices_(vI) { assert(vertexIndices_.size() == N); }
	cell() { vertexIndices_.resize(N); }
	uint& operator ()(uint i) { return vertexIndices_.at(i); }
	uint operator ()(uint i) const { return vertexIndices_.at(i); }
	bool contains(uint i) const { return vertexIndices_.contains(i); }
	uint find(uint i) const { return vertexIndices_.find(i); }
	
	uint size() const { return N; }
	
	friend std::ostream& operator <<(std::ostream& os, const cell& c)
	{
		os << "[";
		for(uint i = 0; i < N; i++) os << c(i) << (i < N-1 ? ", " : "]");
		return os;
	}
	
};




}	// namespace tma

#endif // TMA_GEOMETRY_H
