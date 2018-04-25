#ifndef TMA_QUADRATURE_H_
#define TMA_QUADRATURE_H_

#include <tma/types.h>
#include <functional>


namespace tma
{
// basefunctions phi
template <class T>
struct basefunctions
{
	static uint dim();
	static real phi(uint i, const T& cell, const point<T::D()>& x);
};


// reference cell integration

class TWBQuadrature
{
protected:
	const real lambda1[10] = {0., 1., 0., 0.2673, 0.67282, 0.06492, 0.67165, 0.0654, 0.26938, 0.338674};
	const real lambda2[10] = {1., 0., 0., 0.67282, 0.26733, 0.67165, 0.0649251, 0.26938, 0.06505, 0.33867};
	const real weight[10] = {0.026271, 0.026271, 0.02741, 0.234838, 0.234841, 0.248025, 0.2480305, 0.2518605, 0.2518661, 0.4505789};
	
	point<2> x_q(const triangle& cell, uint q) const
	{
		real l1 = lambda1[q];
		real l2 = lambda2[q];
		real l3 = 1. - l1 - l2;
		real x = cell.a(0)*l1 + cell.b(0)*l2 + cell.c(0)*l3;
		real y = cell.a(1)*l1 + cell.b(1)*l2 + cell.c(1)*l3;
		return point<2>({x, y});
	}
public:
	real Integral(const triangle& cell, real (*f)(const point<2>&)) const
	{
		real s = 0; point<2> p;
	for(uint i = 0; i < 10; i++) {
		s+= f(x_q(cell, i))*weight[i]; 
	}
	return s*cell.area() / 2.;
	}
};

// affine Mapping
/*
template<class T>
struct affineMapping
{
	const T cell;
	std::function<point<T::D()>, const point<T::D()>& p> T_k;
	std::function<QuadraticMatrix<T::D()>, const point<T::D()>& > Jacobian; 
	affineMapping(const T& cell);
};

// affine Mapping Implementation 1D
template<>
affineMapping<interval>::affineMapping(const interval& c) : cell(c)
{
	Jacobian = [&cell] (const point<1>& p) { QuadraticMatrix<1> m; m(0,0) = (p - cell.a); return m; }
	T_k = [&cell] (const point<1>& p) { return cell.a + p(0)*(cell.b-cell.a) };
}

// affine Mapping Implementation 2D
template<>
affineMapping<triangle>::affineMapping(const triangle& cell) : cell(cell)
{
	// find right angle point
	
	
}
*/

// Implementation of basefunctions in 1D interval
template<>
static uint basefunctions<interval>::dim() { return 2; }

template<>
static real basefunctions<interval>::phi(uint i, const interval& I, const point<1>& x)
{
	assert(i == 1 || i == 2);
	if(i == 1)
		return (I.b(0)-x(0))/(I.b(0)-I.a(0));
	return (x(0)-I.a(0))/(I.b(0)-I.a(0));
}

// Implementation of basefunctions in 2D triangle
template<>
static uint basefunctions<triangle>::dim() { return 3; }

template<>
static real basefunctions<triangle>::phi(uint i, const triangle& T, const point<2>& x)
{
	assert(1<= i && i <= 3);
	real l1, l2, d;
	d = (T.b(1)-T.c(1))*(T.a(0) - T.c(0)) + (T.c(1) - T.b(1))*(T.a(1) - T.c(1));
	if(i == 1 || i == 3)
		l1 = ((T.b(1) - T.c(1))*(x(0) - T.c(0)) + (T.c(0)-T.b(0))*(x(1)-T.c(1)))/d;
	if(i == 2 || i == 3)
		l2 = ((T.c(1) - T.a(1))*(x(0) - T.c(0)) + (T.a(0)-T.c(0))*(x(1)-T.c(1)))/d;
	if( i == 1)
		return l1;
	if( i == 2)
		return l2;
	if(i == 3)
		return 1.-l1-l2;
}

// Numerical integration on reference cell


// 

} /* namespace tma */

#endif /* TMA_QUADRATURE_H_ */

