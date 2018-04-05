#ifndef TMA_INTERVAL_H_
#define TMA_INTERVAL_H_

#include <tma/types.h>

namespace tma
{

/*
 * This is a simple definition of an interval cell with reference coordinates.
 */
struct interval
{
  // Return the topological dimension of the cell
  static uint dim() { return 1; }

  // Number of entities of given topological dimension
  static uint num(uint d)
  {
    static uint const _s[2] = { 2, 1 };
    return _s[d];
  }

  struct reference
  {
    // Coordinates of vertices in the reference interval
    static real const * x(uint i)
    {
      static real const _s[2][1] = { { 0.0 }, { 1.0 } };
      return _s[i];
    }
  };
};

} /* namespace tma */

#endif /* TMA_INTERVAL_H_ */
