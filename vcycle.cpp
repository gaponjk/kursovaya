#include "vcycle.h"
#include "gauss_seidel.h"
#include "residual.h"
#include "restriction.h"
#include "prolongation.h"
#include "coefficients.h"

#include <cmath>

void VCycleSolver::init(int N1, int N2,
                        const ProblemParams& params,
                        int numLevels,
                        int nu1, int nu2)
{
    m_numLevels = numLevels;
    m_nu1 = nu1;
    m_nu2 = nu2;
    m_params = params;

    levels.resize(numLevels);

    double h1 = params.l1 / N1;
    double h2 = params.l2 / N2;
    levels[0].allocate(N1, N2, h1, h2);

    for (int lev = 1; lev < numLevels; lev++) {
        int N1_prev = levels[lev - 1].N1;
        int N2_prev = levels[lev - 1].N2;

        int N1_c = N1_prev / 2;
        int N2_c = N2_prev / 2;

        double h1_c = params.l1 / N1_c;
        double h2_c = params.l2 / N2_c;

        levels[lev].allocate(N1_c, N2_c, h1_c, h2_c);
    }
}

void VCycleSolver::prepareCoefficients(double tau, double tj, double tj_prev)
{
    m_tau = tau;
    m_tj = tj;

    for (int lev = 0; lev < m_numLevels; lev++) {
        computeCoefficients(levels[lev], m_params, tau, tj, tj_prev);
    }
}
void VCycleSolver::doVCycle()
{
    for (int lev = 0; lev < m_numLevels - 1; lev++) {
        gaussSeidel(levels[lev], m_nu1);
        computeResidual(levels[lev]);
        restriction(levels[lev], levels[lev + 1]);
        levels[lev + 1].zeroY();
    }

    gaussSeidel(levels[m_numLevels - 1], 200);

    for (int lev = m_numLevels - 2; lev >= 0; lev--) {
        prolongation(levels[lev + 1], levels[lev]);
        gaussSeidel(levels[lev], m_nu2);
    }
}

int VCycleSolver::solve(double tol, int maxCycles)
{
    int cycle;
    for (cycle = 0; cycle < maxCycles; cycle++) {
        doVCycle();
        computeResidual(levels[0]);
        double res = residualNorm(levels[0]);
        if (res < tol) {
            cycle++;  // учтём этот цикл
            break;
        }
    }
    return cycle;
}