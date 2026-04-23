#include "test_problems.h"

double u_exact(double x1, double x2, double t, int testProblem)
{
    if (testProblem == 1) {
        // u = x1^2 + x2^2 + t
        return x1 * x1 + x2 * x2 + t;
    } else {
        // u = exp(x1 + x2 + t)
        return exp(x1 + x2 + t);
    }
}

double k1_func(double x1, double x2, double t, int testProblem)
{
    if (testProblem == 1) {
        (void)x1; (void)x2; (void)t;
        return 1.0;
    } else {
        return exp(x1 + x2 + t);
    }
}

double k2_func(double x1, double x2, double t, int testProblem)
{
    if (testProblem == 1) {
        (void)x1; (void)x2; (void)t;
        return 1.0;
    } else {
        return exp(x1 + x2 + t);
    }
}

double f_rhs(double x1, double x2, double t, int testProblem)
{
    if (testProblem == 1) {
        (void)x1; (void)x2; (void)t;
        return -3.0;
    } else {
        double s = x1 + x2 + t;
        return exp(s) - 4.0 * exp(2.0 * s);
    }
}

double phi(double x1, double x2, int testProblem)
{
    return u_exact(x1, x2, 0.0, testProblem);
}

// x1 = 0:  u(0, x2, t)
double a_left(double x2, double t, int testProblem)
{
    if (testProblem == 1) {
        return x2 * x2 + t;
    } else {
        return exp(x2 + t);
    }
}

double a_right(double x2, double t, double l1, int testProblem)
{
    if (testProblem == 1) {
        return l1 * l1 + x2 * x2 + t;
    } else {
        return exp(l1 + x2 + t);
    }
}

double b_left(double x1, double t, int testProblem)
{
    if (testProblem == 1) {
        return x1 * x1 + t;
    } else {
        return exp(x1 + t);
    }
}

double b_right(double x1, double t, double l2, int testProblem)
{
    if (testProblem == 1) {
        return x1 * x1 + l2 * l2 + t;
    } else {
        return exp(x1 + l2 + t);
    }
}
