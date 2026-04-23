#ifndef TEST_PROBLEMS_H
#define TEST_PROBLEMS_H

#include <cmath>

// Параметры области и задачи
struct ProblemParams
{
    double l1 = 1.0; // длина по x1
    double l2 = 1.0; // длина по x2
    double T  = 1.0; // конечное время
    int testProblem = 1; // 1 или 2
};

// Точное решение u(x1, x2, t)
double u_exact(double x1, double x2, double t, int testProblem);

double k1_func(double x1, double x2, double t, int testProblem);

double k2_func(double x1, double x2, double t, int testProblem);

// Правая часть f(x1, x2, t)
double f_rhs(double x1, double x2, double t, int testProblem);

// Начальное условие u0(x1, x2) = u(x1, x2, 0)
double phi(double x1, double x2, int testProblem);

// Граничные условия (первого рода)

double a_left(double x2, double t, int testProblem);

double a_right(double x2, double t, double l1, int testProblem);

double b_left(double x1, double t, int testProblem);

double b_right(double x1, double t, double l2, int testProblem);

#endif // TEST_PROBLEMS_H