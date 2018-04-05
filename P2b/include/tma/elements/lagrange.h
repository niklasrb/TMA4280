#ifndef TMA_LAGRANGE_H_
#define TMA_LAGRANGE_H_

#include <tma/types.h>

namespace tma
{

struct P1
{
  struct interval
  {
    /// Return the dimension of the reference element
    uint dim() const
    {
      return 2;
    }

    // Coordinates of degrees of freedom in the reference interval
    real const * x(uint i)
    {
      static real const _s[2][1] = { { 0.0 }, { 1.0 } };
      return _s[i];
    }

    // Evaluate the polynomial basis in the reference element
    void operator()(real const* x, real* v) const
    {
      v[0] = 1.0 - x[0];
      v[1] = x[0];
    }

    // Evaluate the polynomial basis derivatives in the reference element
    void d(real const* x, real** v) const
    {
      v[0][0] = - 1.0;
      v[1][0] = + 1.0;
    }
  };

  struct triangle
  {
    /// Return the dimension of the reference element
    uint dim() const
    {
      return 3;
    }

    // Coordinates of degrees of freedom in the reference triangle
    real const * x(uint i)
    {
      static real const x_[3][2] = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 } };
      return x_[i];
    }

    // Evaluate the polynomial basis at x
    void operator()(real const* x, real* v)
    {
      v[0] = 1.0 - x[0] - x[1];
      v[1] = x[0];
      v[2] = x[1];

    }

    // Evaluate the polynomial basis derivatives at x
    void d(real const* x, real** v)
    {
      v[0][0] = -1.0;
      v[1][0] = -1.0;
      v[0][1] = +1.0;
      v[1][1] =  0.0;
      v[0][2] =  0.0;
      v[1][2] = +1.0;
    }
  };

};

} /* namespace tma */

#endif /* TMA_LAGRANGE_H_ */
