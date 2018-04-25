#ifndef TMA_MAPPING_H_
#define TMA_MAPPING_H_

#include <tma/types.h>

namespace tma
{

template<class T>
struct Jacobian
{
  /// Jacobian
  real J[T::D()][T::D()];

  /// Inverse
  real K[T::D()][T::D()];

  /// Determinant
  real det;

  ///
  void operator()(const T& I);

};

//--- IMPLEMENTATION 1D -------------------------------------------------------

template<>
 void Jacobian<interval>::operator()(const interval& I)
{
  // Compute Jacobian
  J[0][0] = I.b(0) - I.a(0);

  // Compute determinant
  det = J[0][0];

  // Compute inverse
  K[0][0] = 1.0 / det;
};

//--- IMPLEMENTATION 2D -------------------------------------------------------

template<>
 void Jacobian<triangle>::operator()(const point<2> x)
{
  // Compute Jacobian
  J[0][0] = x1(0) - x0(0;
  J[1][0] = x[1][1] - x[0][1];
  J[0][1] = x[2][0] - x[0][0];
  J[1][1] = x[2][1] - x[0][1];

  // Compute determinant
  det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

  // Compute inverse
  K[0][0] = + J[1][1] / det;
  K[1][0] = - J[1][0] / det;
  K[0][1] = - J[0][1] / det;
  K[1][1] = + J[0][0] / det;
};

//-----------------------------------------------------------------------------

} /* namespace tma */

#endif /* TMA_MAPPING_H_ */
