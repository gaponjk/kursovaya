#include "restriction.h"


void restriction(const GridLevel& fine, GridLevel& coarse)
{
    int N1c = coarse.N1;
    int N2c = coarse.N2;

    const auto& rf = fine.r;

    for (int i = 0; i <= N1c; i++)
        for (int k = 0; k <= N2c; k++)
            coarse.f[i][k] = 0.0;

    for (int ic = 1; ic <= N1c - 1; ic++) {
        int i_f = 2 * ic;
        for (int kc = 1; kc <= N2c - 1; kc++) {
            int k_f = 2 * kc;

            double center = 4.0 * rf[i_f][k_f];

            double direct = 2.0 * (rf[i_f - 1][k_f] + rf[i_f + 1][k_f]
                                   + rf[i_f][k_f - 1] + rf[i_f][k_f + 1]);

            double diag = rf[i_f - 1][k_f - 1] + rf[i_f + 1][k_f - 1]
                          + rf[i_f - 1][k_f + 1] + rf[i_f + 1][k_f + 1];

            coarse.f[ic][kc] = (center + direct + diag) / 4.0;
        }
    }
}