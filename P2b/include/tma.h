#ifndef TMA_H_
#define TMA_H_

#include <cstdlib>
#include <mpi.h>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <stdint.h>
#include <vector>
#include <map>
#include <iostream>
#include <cassert>
#include <iomanip>
#include <initializer_list>
#include <functional>

#define OPENMP

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
