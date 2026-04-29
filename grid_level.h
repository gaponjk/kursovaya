#ifndef GRID_LEVEL_H
#define GRID_LEVEL_H

#include <vector>


struct GridLevel
{
    int N1;
    int N2;

    double h1;
    double h2;

    std::vector<std::vector<double>> y;

    std::vector<std::vector<double>> cc, ac, bc, cA, cB;

    // правая часть
    std::vector<std::vector<double>> f;

    // невязка
    std::vector<std::vector<double>> r;


    void allocate(int N1_, int N2_, double h1_, double h2_);
    void zeroY();
    void zeroR();
};

#endif // GRID_LEVEL_H