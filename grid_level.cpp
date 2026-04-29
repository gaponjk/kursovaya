#include "grid_level.h"

void GridLevel::allocate(int N1_, int N2_, double h1_, double h2_)
{
    N1 = N1_;
    N2 = N2_;
    h1 = h1_;
    h2 = h2_;

    int rows = N1 + 1;
    int cols = N2 + 1;

    auto make2D = [&]() -> std::vector<std::vector<double>> {
        return std::vector<std::vector<double>>(rows, std::vector<double>(cols, 0.0));
    };

    y = make2D();
    cc = make2D();
    ac = make2D();
    bc = make2D();
    cA = make2D();
    cB = make2D();
    f = make2D();
    r = make2D();
}

void GridLevel::zeroY()
{
    for (auto& row : y)
        std::fill(row.begin(), row.end(), 0.0);
}

void GridLevel::zeroR()
{
    for (auto& row : r)
        std::fill(row.begin(), row.end(), 0.0);
}