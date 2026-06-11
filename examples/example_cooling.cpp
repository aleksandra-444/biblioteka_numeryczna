/**
 * @file example_cooling.cpp
 * @brief Example 1: Solar surface cooling simulation
 *
 * Demonstrates ODE solvers + integration error analysis.
 * Based on metody_numeryczne_9.cpp / metody_numeryczne_10.cpp
 *
 * Physical model: Stefan-Boltzmann cooling law (beta=0, pure radiation)
 *   dT/dt = -alpha * T^4
 *
 * Exact solution:
 *   T(t) = T0*10000 / (1851232043943*t + 1e12)^(1/3)
 *
 * Build:
 *   g++ -std=c++17 -I../include example_cooling.cpp \
 *       ../src/ode.cpp -o example_cooling
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../include/ode.h"

using namespace std;

// ── Physical parameters ───────────────────────────────────────────────────
static const double ALPHA = 3e-12;
static const double T0 = 5903.0;   // K (Sun surface)
static const double T_END = 5903.0;   // seconds simulated
static const int    N_STEPS = 100000;

// ── RHS: dT/dt = -alpha * T^4 ─────────────────────────────────────────────
static RHS f_cool = [](double T, double) { return -ALPHA * T * T * T * T; };

// ── Exact solution ─────────────────────────────────────────────────────────
static double T_exact(double t)
{
    return T0 * 10000.0 / pow(1851232043943.0 * t + 1e12, 1.0 / 3.0);
}

int main()
{
    cout << fixed << setprecision(6);
    cout << "===========================================\n";
    cout << "  EXAMPLE 1: Solar Surface Cooling (ODE)\n";
    cout << "===========================================\n\n";
    cout << "Parameters: T0=" << T0 << " K  |  alpha=" << ALPHA
        << "  |  t_end=" << T_END << " s  |  N=" << N_STEPS << "\n\n";

    // ── Solve with all 4 methods ──────────────────────────────────────────
    double T_e = euler_solve(f_cool, T0, T_END, N_STEPS);
    double T_h = heun_solve(f_cool, T0, T_END, N_STEPS);
    double T_m = midpoint_solve(f_cool, T0, T_END, N_STEPS);
    double T_r = rk4_solve(f_cool, T0, T_END, N_STEPS);
    double T_ex = T_exact(T_END);

    cout << "Final temperature T(t_end):\n";
    cout << setw(20) << "Exact:" << setw(14) << T_ex << " K\n";
    cout << setw(20) << "Euler:" << setw(14) << T_e << " K  |error| = " << fabs(T_e - T_ex) << "\n";
    cout << setw(20) << "Heun:" << setw(14) << T_h << " K  |error| = " << fabs(T_h - T_ex) << "\n";
    cout << setw(20) << "Midpoint:" << setw(14) << T_m << " K  |error| = " << fabs(T_m - T_ex) << "\n";
    cout << setw(20) << "RK4:" << setw(14) << T_r << " K  |error| = " << fabs(T_r - T_ex) << "\n\n";

    // ── MSE over whole trajectory ─────────────────────────────────────────
    function<double(double)> ex_fn = T_exact;
    double mseE = ode_mse(f_cool, ex_fn, T0, T_END, N_STEPS, 0);
    double mseH = ode_mse(f_cool, ex_fn, T0, T_END, N_STEPS, 1);
    double mseM = ode_mse(f_cool, ex_fn, T0, T_END, N_STEPS, 2);
    double mseR = ode_mse(f_cool, ex_fn, T0, T_END, N_STEPS, 3);

    cout << "Mean-Squared Error (MSE) over trajectory:\n";
    cout << setw(20) << "Euler MSE:" << setw(18) << setprecision(10) << mseE << "\n";
    cout << setw(20) << "Heun MSE:" << setw(18) << mseH << "\n";
    cout << setw(20) << "Midpoint MSE:" << setw(18) << mseM << "\n";
    cout << setw(20) << "RK4 MSE:" << setw(18) << mseR << "\n\n";

    // ── Save CSV for plotting ─────────────────────────────────────────────
    ofstream csv("cooling_results.csv");
    csv << "t,Exact,Euler,Heun,Midpoint,RK4\n";
    double h = T_END / N_STEPS;
    double tE = T0, tH = T0, tM = T0, tR = T0, t = 0;
    int print_every = N_STEPS / 50;
    for (int i = 0; i <= N_STEPS; ++i) {
        if (i % print_every == 0 || i == N_STEPS)
            csv << t << "," << T_exact(t) << ","
            << tE << "," << tH << "," << tM << "," << tR << "\n";
        if (i < N_STEPS) {
            tE = euler_step(f_cool, tE, t, h);
            tH = heun_step(f_cool, tH, t, h);
            tM = midpoint_step(f_cool, tM, t, h);
            tR = rk4_step(f_cool, tR, t, h);
            t += h;
        }
    }
    csv.close();
    cout << "Results saved to: cooling_results.csv\n";

    return 0;
}