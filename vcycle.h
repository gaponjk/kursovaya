#ifndef VCYCLE_H
#define VCYCLE_H

#include "grid_level.h"
#include "test_problems.h"
#include <vector>

class VCycleSolver
{
public:
    void init(int N1, int N2,
              const ProblemParams& params,
              int numLevels,
              int nu1, int nu2);

    void prepareCoefficients(double tau, double tj, double tj_prev);

    void doVCycle();

    int solve(double tol, int maxCycles);
    GridLevel& finest() { return levels[0]; }
    const GridLevel& finest() const { return levels[0]; }

    std::vector<GridLevel> levels;

private:
    int m_numLevels;
    int m_nu1;
    int m_nu2;
    ProblemParams m_params;
    double m_tau;
    double m_tj;
};

#endif // VCYCLE_H