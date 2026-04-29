#include "coefficients.h"
#include <cmath>

void computeCoefficients(GridLevel& level, const ProblemParams& params,
                         double tau, double tj, double tj_prev)
{
    (void)tj_prev;
    int N1 = level.N1;
    int N2 = level.N2;
    double h1 = level.h1;
    double h2 = level.h2;

    for (int i = 1; i <= N1 - 1; i++) {
        double x1 = i * h1;
        for (int k = 1; k <= N2 - 1; k++) {
            double x2 = k * h2;

            double k1_plus  = k1_func(x1 + 0.5*h1, x2, tj, params.testProblem);
            double k1_minus = k1_func(x1 - 0.5*h1, x2, tj, params.testProblem);

            double k2_plus  = k2_func(x1, x2 + 0.5*h2, tj, params.testProblem);
            double k2_minus = k2_func(x1, x2 - 0.5*h2, tj, params.testProblem);

            level.ac[i][k] = -k1_minus;
            level.bc[i][k] = -k1_plus;


            level.cA[i][k] = -(h1*h1/(h2*h2)) * k2_minus;
            level.cB[i][k] = -(h1*h1/(h2*h2)) * k2_plus;
            level.cc[i][k] = h1*h1/tau
                             + k1_plus + k1_minus
                             + (h1*h1/(h2*h2)) * (k2_plus + k2_minus);
        }
    }
}


void computeRHS(GridLevel& level, const ProblemParams& params,
                double tau, double tj, double tj_prev,
                const std::vector<std::vector<double>>& y_prev)
{
    (void)tj_prev;
    int N1 = level.N1;
    int N2 = level.N2;
    double h1 = level.h1;
    double h2 = level.h2;

    for (int k = 0; k <= N2; k++) {
        double x2 = k * h2;
        level.y[0][k]  = u_exact(0.0,       x2, tj, params.testProblem);
        level.y[N1][k] = u_exact(params.l1, x2, tj, params.testProblem);
    }
    for (int i = 0; i <= N1; i++) {
        double x1 = i * h1;
        level.y[i][0]  = u_exact(x1, 0.0,        tj, params.testProblem);
        level.y[i][N2] = u_exact(x1, params.l2, tj, params.testProblem);
    }
    for (int i = 1; i <= N1 - 1; i++) {
        double x1 = i * h1;
        for (int k = 1; k <= N2 - 1; k++) {
            double x2 = k * h2;
            level.f[i][k] = h1 * h1 / tau * y_prev[i][k]
                            + h1 * h1 * f_rhs(x1, x2, tj, params.testProblem);
        }
    }
}
