#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include "grid_level.h"
#include "test_problems.h"

void computeCoefficients(GridLevel& level,
                         const ProblemParams& params,
                         double tau,
                         double tj,
                         double tj_prev
                         );

void computeRHS(GridLevel& level,
                const ProblemParams& params,
                double tau,
                double tj,
                double tj_prev,
                const std::vector<std::vector<double>>& y_prev
                );

#endif // COEFFICIENTS_H