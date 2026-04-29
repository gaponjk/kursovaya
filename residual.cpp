#include "residual.h"
#include <cmath>
#include <algorithm>

void computeResidual(GridLevel& level)
{
    int N1 = level.N1;
    int N2 = level.N2;

    auto& y  = level.y;
    auto& cc = level.cc;
    auto& ac = level.ac;
    auto& bc = level.bc;
    auto& cA = level.cA;
    auto& cB = level.cB;
    auto& f  = level.f;
    auto& r  = level.r;

    for (int k = 0; k <= N2; k++) {
        r[0][k]  = 0.0;
        r[N1][k] = 0.0;
    }
    for (int i = 0; i <= N1; i++) {
        r[i][0]  = 0.0;
        r[i][N2] = 0.0;
    }

    for (int i = 1; i <= N1 - 1; i++) {
        for (int k = 1; k <= N2 - 1; k++) {
            double Ay = cc[i][k] * y[i][k]
                        + ac[i][k] * y[i][k-1]
                        + bc[i][k] * y[i][k+1]
                        + cA[i][k] * y[i-1][k]
                        + cB[i][k] * y[i+1][k];
            r[i][k] = f[i][k] - Ay;
        }
    }
}

double residualNorm(const GridLevel& level)
{
    int N1 = level.N1;
    int N2 = level.N2;
    double maxR = 0.0;

    for (int i = 1; i <= N1 - 1; i++) {
        for (int k = 1; k <= N2 - 1; k++) {
            maxR = std::max(maxR, std::fabs(level.r[i][k]));
        }
    }
    return maxR;
}