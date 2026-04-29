#include "gauss_seidel.h"

void gaussSeidel(GridLevel& level, int numIterations)
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

    for (int iter = 0; iter < numIterations; iter++) {
        for (int i = 1; i <= N1 - 1; i++) {
            for (int k = 1; k <= N2 - 1; k++) {
                y[i][k] = ( f[i][k]
                           - ac[i][k] * y[i][k-1]
                           - bc[i][k] * y[i][k+1]
                           - cA[i][k] * y[i-1][k]
                           - cB[i][k] * y[i+1][k]
                           ) / cc[i][k];
            }
        }
    }
}