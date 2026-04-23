#include "prolongation.h"

void prolongation(const GridLevel& coarse, GridLevel& fine)
{
    int N1c = coarse.N1;
    int N2c = coarse.N2;
    int N1f = fine.N1;
    int N2f = fine.N2;

    const auto& ec = coarse.y;


    for (int ic = 1; ic <= N1c - 1; ic++) {
        int i_f = 2 * ic;
        for (int kc = 1; kc <= N2c - 1; kc++) {
            int k_f = 2 * kc;
            fine.y[i_f][k_f] += ec[ic][kc];
        }
    }

    for (int ic = 0; ic <= N1c - 1; ic++) {
        int i_f = 2 * ic + 1;
        if (i_f < 1 || i_f > N1f - 1) continue;
        for (int kc = 1; kc <= N2c - 1; kc++) {
            int k_f = 2 * kc;
            fine.y[i_f][k_f] += 0.5 * (ec[ic][kc] + ec[ic + 1][kc]);
        }
    }

    for (int ic = 1; ic <= N1c - 1; ic++) {
        int i_f = 2 * ic;
        for (int kc = 0; kc <= N2c - 1; kc++) {
            int k_f = 2 * kc + 1;
            if (k_f < 1 || k_f > N2f - 1) continue;
            fine.y[i_f][k_f] += 0.5 * (ec[ic][kc] + ec[ic][kc + 1]);
        }
    }

    for (int ic = 0; ic <= N1c - 1; ic++) {
        int i_f = 2 * ic + 1;
        if (i_f < 1 || i_f > N1f - 1) continue;
        for (int kc = 0; kc <= N2c - 1; kc++) {
            int k_f = 2 * kc + 1;
            if (k_f < 1 || k_f > N2f - 1) continue;
            fine.y[i_f][k_f] += 0.25 * (ec[ic][kc]   + ec[ic + 1][kc]
                                        + ec[ic][kc+1] + ec[ic + 1][kc + 1]);
        }
    }
}