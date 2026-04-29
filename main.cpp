#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>

#include "grid_level.h"
#include "coefficients.h"
#include "gauss_seidel.h"
#include "residual.h"
#include "vcycle.h"
#include "test_problems.h"

// Вычислить ошибку (максимум-норма) относительно точного решения
double computeError(const GridLevel& level, double t,
                    const ProblemParams& params)
{
    double maxErr = 0.0;
    for (int i = 0; i <= level.N1; i++) {
        double x1 = i * level.h1;
        for (int k = 0; k <= level.N2; k++) {
            double x2 = k * level.h2;
            double exact = u_exact(x1, x2, t, params.testProblem);
            double err = std::fabs(level.y[i][k] - exact);
            if (err > maxErr) maxErr = err;
        }
    }
    return maxErr;
}

// Задание начального условия на мелкой сетке
void setInitialCondition(GridLevel& level, const ProblemParams& params)
{
    for (int i = 0; i <= level.N1; i++) {
        double x1 = i * level.h1;
        for (int k = 0; k <= level.N2; k++) {
            double x2 = k * level.h2;
            level.y[i][k] = phi(x1, x2, params.testProblem);
        }
    }
}
// Решение с помощью чистого Гаусса-Зейделя (для сравнения)
void solveGaussSeidel(int N1, int N2, const ProblemParams& params,
                      double tau, int Nt)
{
    std::cout << "========================================" << std::endl;
    std::cout << "Gauss-Seidel solver" << std::endl;
    std::cout << "N1=" << N1 << " N2=" << N2
              << " tau=" << tau << " Nt=" << Nt << std::endl;
    std::cout << "Test problem: " << params.testProblem << std::endl;
    std::cout << "========================================" << std::endl;

    double h1 = params.l1 / N1;
    double h2 = params.l2 / N2;

    GridLevel level;
    level.allocate(N1, N2, h1, h2);

    setInitialCondition(level, params);

    auto y_prev = level.y;

    auto t_start = std::chrono::high_resolution_clock::now();

    long long totalGSIters = 0;

    for (int j = 1; j <= Nt; j++) {
        double tj      = j * tau;
        double tj_prev = (j - 1) * tau;

        computeCoefficients(level, params, tau, tj, tj_prev);
        computeRHS(level, params, tau, tj, tj_prev, y_prev);

        double tol = 1e-10;
        int maxIter = 20000;
        int iters = 0;

        for (int iter = 1; iter <= maxIter; iter++) {
            gaussSeidel(level, 1);
            computeResidual(level);
            double norm = residualNorm(level);
            iters = iter;
            if (norm < tol) break;
        }

        totalGSIters += iters;

        if (j == 1 || j == Nt || j % std::max(1, Nt / 5) == 0) {
            computeResidual(level);
            double norm = residualNorm(level);
            double error = computeError(level, tj, params);
            std::cout << "  t=" << std::fixed << std::setprecision(4) << tj
                      << "  GS iters=" << std::setw(6) << iters
                      << "  ||r||=" << std::scientific << std::setprecision(2) << norm
                      << "  error=" << error
                      << std::endl;
        }

        y_prev = level.y;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    double error = computeError(level, Nt * tau, params);
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Final error     = " << std::scientific
              << std::setprecision(6) << error << std::endl;
    std::cout << "Total GS iters  = " << totalGSIters << std::endl;
    std::cout << "Time            = " << std::fixed
              << std::setprecision(3) << elapsed << " s" << std::endl;
    std::cout << std::endl;
}

// Решение с помощью многосеточного V-цикла
void solveMultigrid(int N1, int N2, const ProblemParams& params,
                    double tau, int Nt,
                    int numLevels, int nu1, int nu2)
{
    std::cout << "========================================" << std::endl;
    std::cout << "Multigrid V-cycle solver" << std::endl;
    std::cout << "N1=" << N1 << " N2=" << N2
              << " tau=" << tau << " Nt=" << Nt << std::endl;
    std::cout << "Levels=" << numLevels
              << " nu1=" << nu1 << " nu2=" << nu2 << std::endl;
    std::cout << "Test problem: " << params.testProblem << std::endl;
    std::cout << "========================================" << std::endl;

    // Печатаем иерархию сеток
    {
        int n1 = N1, n2 = N2;
        for (int lev = 0; lev < numLevels; lev++) {
            int unknowns = (n1 - 1) * (n2 - 1);
            std::cout << "  Level " << lev
                      << ": N1=" << n1 << " N2=" << n2
                      << "  unknowns=" << unknowns << std::endl;
            n1 /= 2;
            n2 /= 2;
        }
    }

    VCycleSolver mg;
    mg.init(N1, N2, params, numLevels, nu1, nu2);

    GridLevel& finest = mg.finest();
    setInitialCondition(finest, params);

    auto y_prev = finest.y;

    auto t_start = std::chrono::high_resolution_clock::now();

    int totalCycles = 0;

    for (int j = 1; j <= Nt; j++) {
        double tj      = j * tau;
        double tj_prev = (j - 1) * tau;

        mg.prepareCoefficients(tau, tj, tj_prev);
        computeRHS(finest, params, tau, tj, tj_prev, y_prev);

        double tol = 1e-10;
        int maxCycles = 50;
        int cycles = mg.solve(tol, maxCycles);
        totalCycles += cycles;

        if (j == 1 || j == Nt || j % std::max(1, Nt / 5) == 0) {
            double error = computeError(finest, tj, params);
            computeResidual(finest);
            double norm = residualNorm(finest);
            std::cout << "  t=" << std::fixed << std::setprecision(4) << tj
                      << "  V-cycles=" << std::setw(3) << cycles
                      << "  ||r||=" << std::scientific << std::setprecision(2) << norm
                      << "  error=" << error
                      << std::endl;
        }

        y_prev = finest.y;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    double error = computeError(finest, Nt * tau, params);
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Final error     = " << std::scientific
              << std::setprecision(6) << error << std::endl;
    std::cout << "Total V-cycles  = " << totalCycles << std::endl;
    std::cout << "Time            = " << std::fixed
              << std::setprecision(3) << elapsed << " s" << std::endl;
    std::cout << std::endl;
}

// Сравнительная таблица сходимости: MG vs GS на разных сетках
void convergenceStudy(const ProblemParams& params,
                      int numLevelsRequested, int nu1, int nu2,double tau)
{
    std::cout << "============================================" << std::endl;
    std::cout << "Convergence study: Test problem "
              << params.testProblem << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << std::setw(6)  << "N"
              << std::setw(14) << "Error_MG"
              << std::setw(14) << "Error_GS"
              << std::setw(10) << "MGcycles"
              << std::setw(10) << "GSiters"
              << std::setw(12) << "Time_MG"
              << std::setw(12) << "Time_GS"
              << std::endl;
    std::cout << std::string(78, '-') << std::endl;

    std::vector<int> sizes = {16, 32, 64, 128};

    for (int N : sizes) {
        int N1 = N;
        int N2 = N;
        int Nt = static_cast<int>(params.T / tau);
        if (Nt < 1) Nt = 1;

        int maxLev = numLevelsRequested;
        {
            int n = N;
            for (int l = 1; l < maxLev; l++) {
                if (n < 4 || n % 2 != 0) {
                    maxLev = l;
                    break;
                }
                n /= 2;
            }
        }

        double h1 = params.l1 / N1;
        double h2 = params.l2 / N2;

        VCycleSolver mg;
        mg.init(N1, N2, params, maxLev, nu1, nu2);

        GridLevel& finest_mg = mg.finest();
        setInitialCondition(finest_mg, params);
        auto y_prev_mg = finest_mg.y;

        auto t1 = std::chrono::high_resolution_clock::now();
        int totalCycles = 0;
        for (int j = 1; j <= Nt; j++) {
            double tj      = j * tau;
            double tj_prev = (j - 1) * tau;
            mg.prepareCoefficients(tau, tj, tj_prev);
            computeRHS(finest_mg, params, tau, tj, tj_prev, y_prev_mg);
            totalCycles += mg.solve(1e-10, 50);
            y_prev_mg = finest_mg.y;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        double timeMG  = std::chrono::duration<double>(t2 - t1).count();
        double errorMG = computeError(finest_mg, Nt * tau, params);

        GridLevel level_gs;
        level_gs.allocate(N1, N2, h1, h2);
        setInitialCondition(level_gs, params);
        auto y_prev_gs = level_gs.y;

        t1 = std::chrono::high_resolution_clock::now();
        long long totalGS = 0;
        for (int j = 1; j <= Nt; j++) {
            double tj      = j * tau;
            double tj_prev = (j - 1) * tau;
            computeCoefficients(level_gs, params, tau, tj, tj_prev);
            computeRHS(level_gs, params, tau, tj, tj_prev, y_prev_gs);
            for (int iter = 1; iter <= 20000; iter++) {
                gaussSeidel(level_gs, 1);
                computeResidual(level_gs);
                double norm = residualNorm(level_gs);
                totalGS++;
                if (norm < 1e-10) break;
            }
            y_prev_gs = level_gs.y;
        }
        t2 = std::chrono::high_resolution_clock::now();
        double timeGS  = std::chrono::duration<double>(t2 - t1).count();
        double errorGS = computeError(level_gs, Nt * tau, params);

        std::cout << std::setw(6)  << N
                  << std::setw(14) << std::scientific << std::setprecision(4) << errorMG
                  << std::setw(14) << errorGS
                  << std::setw(10) << totalCycles
                  << std::setw(10) << totalGS
                  << std::setw(12) << std::fixed << std::setprecision(3) << timeMG
                  << std::setw(12) << timeGS
                  << std::endl;
    }
    std::cout << std::endl;
}

void levelsStudy(const ProblemParams& params, int N, double tau)
{
    int Nt = static_cast<int>(params.T / tau);
    if (Nt < 1) Nt = 1;

    std::cout << "============================================" << std::endl;
    std::cout << "Effect of number of levels (Problem "
              << params.testProblem << ", N=" << N << ")" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << std::setw(8)  << "Levels"
              << std::setw(14) << "Error"
              << std::setw(12) << "Cycles"
              << std::setw(12) << "Time"
              << std::endl;
    std::cout << std::string(46, '-') << std::endl;

    for (int numLev = 2; numLev <= 6; numLev++) {
        int n = N;
        bool ok = true;
        for (int l = 1; l < numLev; l++) {
            if (n < 4 || n % 2 != 0) { ok = false; break; }
            n /= 2;
        }
        if (!ok) continue;

        double h1 = params.l1 / N;
        double h2 = params.l2 / N;
        (void)h1; (void)h2;

        VCycleSolver mg;
        mg.init(N, N, params, numLev, 2, 2);

        GridLevel& finest = mg.finest();
        setInitialCondition(finest, params);
        auto y_prev = finest.y;

        auto t1 = std::chrono::high_resolution_clock::now();
        int totalCycles = 0;
        for (int j = 1; j <= Nt; j++) {
            double tj      = j * tau;
            double tj_prev = (j - 1) * tau;
            mg.prepareCoefficients(tau, tj, tj_prev);
            computeRHS(finest, params, tau, tj, tj_prev, y_prev);
            totalCycles += mg.solve(1e-10, 50);
            y_prev = finest.y;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        double timeMG = std::chrono::duration<double>(t2 - t1).count();
        double error  = computeError(finest, Nt * tau, params);

        std::cout << std::setw(8)  << numLev
                  << std::setw(14) << std::scientific << std::setprecision(4) << error
                  << std::setw(12) << totalCycles
                  << std::setw(12) << std::fixed << std::setprecision(3) << timeMG
                  << std::endl;
    }
    std::cout << std::endl;
}

void smoothingStudy(const ProblemParams& params, int N, double tau,
                    int numLevels)
{
    int Nt = static_cast<int>(params.T / tau);
    if (Nt < 1) Nt = 1;

    std::cout << "============================================" << std::endl;
    std::cout << "Effect of smoothing steps (Problem "
              << params.testProblem << ", N=" << N << ")" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << std::setw(6)  << "nu1"
              << std::setw(6)  << "nu2"
              << std::setw(14) << "Error"
              << std::setw(12) << "Cycles"
              << std::setw(12) << "Time"
              << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    std::vector<std::pair<int,int>> smoothings = {
        {1,1}, {2,2}, {3,3}, {4,4}, {2,1}, {1,2}
    };

    for (auto& s : smoothings) {
        int nu1 = s.first;
        int nu2 = s.second;

        VCycleSolver mg;
        mg.init(N, N, params, numLevels, nu1, nu2);

        GridLevel& finest = mg.finest();
        setInitialCondition(finest, params);
        auto y_prev = finest.y;

        auto t1 = std::chrono::high_resolution_clock::now();
        int totalCycles = 0;
        for (int j = 1; j <= Nt; j++) {
            double tj      = j * tau;
            double tj_prev = (j - 1) * tau;
            mg.prepareCoefficients(tau, tj, tj_prev);
            computeRHS(finest, params, tau, tj, tj_prev, y_prev);
            totalCycles += mg.solve(1e-10, 50);
            y_prev = finest.y;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        double timeMG = std::chrono::duration<double>(t2 - t1).count();
        double error  = computeError(finest, Nt * tau, params);

        std::cout << std::setw(6)  << nu1
                  << std::setw(6)  << nu2
                  << std::setw(14) << std::scientific << std::setprecision(4) << error
                  << std::setw(12) << totalCycles
                  << std::setw(12) << std::fixed << std::setprecision(3) << timeMG
                  << std::endl;
    }
    std::cout << std::endl;
}

int main()
{
    std::cout << std::setprecision(6);

    //  ТЕСТОВАЯ ЗАДАЧА 1
    {
        ProblemParams params;
        params.l1 = 1.0;
        params.l2 = 1.0;
        params.T  = 0.1;
        params.testProblem = 1;

        int N1 = 64;
        int N2 = 64;
        double tau = 0.01;
        int Nt = static_cast<int>(params.T / tau);

        int numLevels = 4; // 64 -> 32 -> 16 -> 8
        int nu1 = 2;
        int nu2 = 2;

        std::cout << "============================================" << std::endl;
        std::cout << "  TEST PROBLEM 1                            " << std::endl;
        std::cout << "  k1 = k2 = 1                               " << std::endl;
        std::cout << "  u(x1,x2,t) = x1^2 + x2^2 + t              " << std::endl;
        std::cout << "  f = -3                                    " << std::endl;
        std::cout << "============================================" << std::endl;
        std::cout << std::endl;

        // Подробный запуск MG
        solveMultigrid(N1, N2, params, tau, Nt, numLevels, nu1, nu2);

        // Подробный запуск GS
        solveGaussSeidel(N1, N2, params, tau, Nt);

        // Таблица сходимости на разных сетках
        convergenceStudy(params, numLevels, nu1, nu2, tau);

        // Влияние числа уровней
        levelsStudy(params, N1, tau);

        // Влияние числа сглаживаний
        smoothingStudy(params, N1, tau, numLevels);
    }

    //  ТЕСТОВАЯ ЗАДАЧА 2
    {
        ProblemParams params;
        params.l1 = 1.0;
        params.l2 = 1.0;
        params.T  = 0.1;
        params.testProblem = 2;

        int N1 = 64;
        int N2 = 64;
        double tau = 0.01;
        int Nt = static_cast<int>(params.T / tau);

        int numLevels = 4;
        int nu1 = 2;
        int nu2 = 2;

        std::cout << "============================================" << std::endl;
        std::cout << "  TEST PROBLEM 2                            " << std::endl;
        std::cout << "  k1 = k2 = exp(x1+x2+t)                    " << std::endl;
        std::cout << "  u(x1,x2,t) = exp(x1+x2+t)                 " << std::endl;
        std::cout << "  f = exp(x1+x2+t) - 4*exp(2*(x1+x2+t))     " << std::endl;
        std::cout << "============================================" << std::endl;
        std::cout << std::endl;

        // Подробный запуск MG
        solveMultigrid(N1, N2, params, tau, Nt, numLevels, nu1, nu2);

        // Подробный запуск GS
        solveGaussSeidel(N1, N2, params, tau, Nt);

        // Таблица сходимости
        convergenceStudy(params, numLevels, nu1, nu2, tau);

        // Влияние числа уровней
        levelsStudy(params, N1, tau);

        // Влияние числа сглаживаний
        smoothingStudy(params, N1, tau, numLevels);
    }

    std::cout << "  All tests completed successfully.         " << std::endl;

    return 0;
}