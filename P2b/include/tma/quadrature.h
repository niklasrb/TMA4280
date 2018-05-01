#ifndef TMA_QUADRATURE_H_
#define TMA_QUADRATURE_H_

#include <tma/types.h>



namespace tma
{
// basefunctions phi
template <class T>
struct basefunctions
{
	 uint dim() const;
	 real phi(uint i, const T& cell, const point<T::D()>& x) const;
};


// reference cell integration in 2D

class TWBQuadrature
{
protected:
	const real lambda1[10] = {0., 1., 0., 0.2673273531185, 0.6728175529461, 0.0649236350054, 0.6716498539042, 0.0654032456800, 0.2693767069140, 0.3386738503896};
	const real lambda2[10] = {1., 0., 0., 0.6728199218710, 0.2673288599482, 0.6716530111494, 0.0649251690029, 0.2693789366453, 0.0654054874919, 0.3386799893027};
	const real weight[10] = {0.0262712099504, 0.0262716612068, 0.0274163947600, 0.2348383865823, 0.2348412238268, 0.2480251793114, 0.2480304922521, 0.2518604605529, 0.2518660533658, 0.4505789381914};
	
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
	real Integral(const triangle& cell, const std::function<real(const point<2>&)>& f) const
	{
		real s = 0; point<2> p;
		for(uint i = 0; i < 10; i++) 
			s+= f(x_q(cell, i))*weight[i]; 
		return s*cell.area() / 2.;
	}
};


// reference cell integration in 1D
class GaussianQuadrature
{
constexpr static uint  n = 5;
	real x[n] = {0., 0.538469, -0.538469, 0.90618, -0.90618};
	real weight[n] = {0.568889,  0.478629, 0.478629,  0.236927, 0.236927};
	
public:
	real Integral(const interval& cell, const std::function<real(const point<1>&)>& f) const
	{
		real s = 0; point<1> p;
		for(uint i = 0; i < n; i++)
			s += weight[i] * f( (cell.b(0)-cell.a(0))/2. * x[i] + (cell.a(0)+cell.b(0))/2.);
		
		return s*cell.area()/2.;
	}
	
};


// Implementation of basefunctions in 1D interval
template<>
uint basefunctions<interval>::dim() const { return 2; }

template<>
real basefunctions<interval>::phi(uint i, const interval& I, const point<1>& x) const
{
	assert(i < 2);
	if(x(0) < I.a(0) || x(0) > I.b(0))
		return 0.;
	if(i == 0)
		return (I.b(0)-x(0))/(I.b(0)-I.a(0));
	return (x(0)-I.a(0))/(I.b(0)-I.a(0));
}

// Implementation of basefunctions in 2D triangle
template<>
uint basefunctions<triangle>::dim()  const { return 3; }

template<>
real basefunctions<triangle>::phi(uint i, const triangle& T, const point<2>& x) const
{
	assert(i < 3);
	real l1, l2, l3, d;
	d = (T.b(1)-T.c(1))*(T.a(0) - T.c(0)) + (T.c(0) - T.b(0))*(T.a(1) - T.c(1));
	l1 = ((T.b(1) - T.c(1))*(x(0) - T.c(0)) + (T.c(0)-T.b(0))*(x(1)-T.c(1)))/d;
	l2 = ((T.c(1) - T.a(1))*(x(0) - T.c(0)) + (T.a(0)-T.c(0))*(x(1)-T.c(1)))/d;
	l3 = 1. - l1 - l2;
	if(l1 < 0 || l1 > 1 || l2 < 0 || l2 > 1 || l3 <0 || l3 > 1)	// outside triangle
		return 0.;
	if( i == 0)
		return l1;
	if( i == 1)
		return l2;
	return l3;
}

} /* namespace tma */

#endif /* TMA_QUADRATURE_H_ */

