#ifndef RESIDUAL_H
#define RESIDUAL_H

#include "grid_level.h"


void computeResidual(GridLevel& level);

double residualNorm(const GridLevel& level);

#endif // RESIDUAL_H