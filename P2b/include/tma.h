#ifndef TMA_H_
#define TMA_H_

#include "tma/test.h"

#include <cstdlib>
#include <mpi.h>
#include <cmath>

#ifdef OPENMP
#include <omp.h>
#endif

#include "tma/geometry.h"
#include "tma/mesh/mesh.h"

#include "tma/data/matrix.h"
#include "tma/data/vector.h"

#include "tma/quadrature.h"
#include "tma/helmholtz.h"
#include "tma/parallel.h"
#include "tma/uheat.h"


#endif /* TMA_H_ */
