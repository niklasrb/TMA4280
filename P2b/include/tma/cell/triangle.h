#ifndef TMA_TRIANGLE_H_
#define TMA_TRIANGLE_H_

#include <tma/types.h>

namespace tma
{

/*
 * This is a simple definition of a triangle cell with reference coordinates.
 */
struct triangle
{
  // Return the topological dimension of the cell
  static uint dim() { return 2; }

  // Number of entities of given topological dimension
  static uint num(uint d)
  {
    static uint const _s[3] = { 3, 3, 1 };
    return _s[d];
  }

  struct reference
  {
    // Coordinates of vertices in the reference triangle
    static real const * x(uint i)
    {
      static real const x_[3][2] = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 } };
      return x_[i];
    }
  };
};

} /* namespace tma */

#endif /* TMA_TRIANGLE_H_ */
