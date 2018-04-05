#ifndef TMA_MAPPING_H_
#define TMA_MAPPING_H_

#include <tma/types.h>

namespace tma
{

template<uint N>
struct Jacobian
{
  /// Jacobian
  real J[N][N];

  /// Inverse
  real K[N][N];

  /// Determinant
  real det;

  ///
  void operator()(real const * const * x);

};

//--- IMPLEMENTATION 1D -------------------------------------------------------

template<>
inline void Jacobian<1>::operator()(real const * const * x)
{
  // Compute Jacobian
  J[0][0] = x[1][0] - x[0][0];

  // Compute determinant
  det = J[0][0];

  // Compute inverse
  K[0][0] = 1.0 / det;
};

//--- IMPLEMENTATION 2D -------------------------------------------------------

template<>
inline void Jacobian<2>::operator()(real const * const * x)
{
  // Compute Jacobian
  J[0][0] = x[1][0] - x[0][0];
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
