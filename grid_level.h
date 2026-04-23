#ifndef GRID_LEVEL_H
#define GRID_LEVEL_H

#include <vector>


struct GridLevel
{
    int N1;
    int N2;
    // внутренние узлы по x2: 1..N2-1

    double h1; // шаг по x1
    double h2; // шаг по x2

    std::vector<std::vector<double>> y;

    std::vector<std::vector<double>> cc, ac, bc, cA, cB;

    // f[i][k] — правая часть
    std::vector<std::vector<double>> f;

    // r[i][k] — невязка
    std::vector<std::vector<double>> r;

    std::vector<std::vector<double>> u_saved;
    std::vector<std::vector<double>> f_saved;

    void allocate(int N1_, int N2_, double h1_, double h2_);
    void zeroY();
    void zeroR();
};

#endif // GRID_LEVEL_H