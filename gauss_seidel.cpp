#include "gauss_seidel.h"

void gaussSeidel(GridLevel& level, int l_max)
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

    for (int l = 0; l < l_max; l++) {
        for (int i = 1; i <= N1 - 1; i++) {
            for (int k = 1; k <= N2 - 1; k++) {
                y[i][k] = (f[i][k]
                           - ac[i][k] * y[i-1][k]
                           - bc[i][k] * y[i+1][k]
                           - cA[i][k] * y[i][k-1]
                           - cB[i][k] * y[i][k+1]) / cc[i][k];
            }
        }
    }
}