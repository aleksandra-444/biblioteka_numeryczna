/**
 * @file test_ode.cpp
 * @brief Unit tests for ODE solvers
 *
 * Build:  g++ -std=c++17 -I../include test_ode.cpp ../src/ode.cpp -o test_ode
 * Run:    ./test_ode
 */

#include <iostream>
#include <cmath>
#include "ode.h"

using namespace std;

static int passed = 0, failed = 0;

static void check(bool cond, const string& name)
{
    if (cond) { cout << "  [PASS] " << name << "\n"; ++passed; }
    else { cout << "  [FAIL] " << name << "\n"; ++failed; }
}

// Exponential decay: dy/dt = -y,  y(0)=1  =>  y(t)=exp(-t)
static auto f_decay = [](double y, double /*t*/) { return -y; };
static auto ex_decay = [](double t) { return exp(-t); };

// Linear growth: dy/dt = 2t,  y(0)=0  =>  y(t)=t^2
static auto f_linear = [](double /*y*/, double t) { return 2.0 * t; };
static auto ex_linear = [](double t) { return t * t; };

// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
void test_euler()
{
    cout << "\n[euler_solve]:\n";
    // y(1) should be near exp(-1) ~ 0.3679
    double y = euler_solve(f_decay, 1.0, 1.0, 10000);
    check(fabs(y - exp(-1.0)) < 1e-4, "decay y(1) near exp(-1), tol 1e-4");

    // Linear RHS: y(3) = 9
    double y2 = euler_solve(f_linear, 0.0, 3.0, 10000);
    check(fabs(y2 - 9.0) < 1e-3, "y'=2t, y(3)=9, tol 1e-3");
}

void test_heun()
{
    cout << "\n[heun_solve]:\n";
    double y = heun_solve(f_decay, 1.0, 1.0, 1000);
    check(fabs(y - exp(-1.0)) < 1e-6, "Heun: decay y(1) near exp(-1), tol 1e-6");

    double y2 = heun_solve(f_linear, 0.0, 3.0, 1000);
    check(fabs(y2 - 9.0) < 1e-6, "Heun: y'=2t, y(3)=9, tol 1e-6");
}

void test_midpoint()
{
    cout << "\n[midpoint_solve]:\n";
    double y = midpoint_solve(f_decay, 1.0, 1.0, 1000);
    check(fabs(y - exp(-1.0)) < 1e-6, "Midpoint: decay y(1), tol 1e-6");

    double y2 = midpoint_solve(f_linear, 0.0, 3.0, 1000);
    check(fabs(y2 - 9.0) < 1e-6, "Midpoint: y'=2t, y(3)=9, tol 1e-6");
}

void test_rk4()
{
    cout << "\n[rk4_solve]:\n";
    // RK4 should be very precise even with few steps
    double y = rk4_solve(f_decay, 1.0, 1.0, 100);
    check(fabs(y - exp(-1.0)) < 1e-10, "RK4: decay y(1), tol 1e-10");

    double y2 = rk4_solve(f_linear, 0.0, 3.0, 100);
    check(fabs(y2 - 9.0) < 1e-10, "RK4: y'=2t, y(3)=9, tol 1e-10");
}

// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
// Accuracy ordering: Euler < Heun ~ Midpoint < RK4
// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
void test_accuracy_ordering()
{
    cout << "\n[accuracy ordering] RK4 > Heun > Euler (same N):\n";
    int N = 50;
    double yE = euler_solve(f_decay, 1.0, 1.0, N);
    double yH = heun_solve(f_decay, 1.0, 1.0, N);
    double yRK = rk4_solve(f_decay, 1.0, 1.0, N);
    double exact = exp(-1.0);

    check(fabs(yRK - exact) < fabs(yH - exact), "RK4 more accurate than Heun");
    check(fabs(yH - exact) < fabs(yE - exact), "Heun more accurate than Euler");
}

// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
// ode_mse: RK4 MSE must be orders of magnitude lower than Euler
// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
void test_ode_mse()
{
    cout << "\n[ode_mse]:\n";
    double mseE = ode_mse(f_decay, ex_decay, 1.0, 1.0, 200, 0);
    double mseRK = ode_mse(f_decay, ex_decay, 1.0, 1.0, 200, 3);
    check(mseRK < mseE, "RK4 MSE < Euler MSE");
    check(mseRK < 1e-15, "RK4 MSE < 1e-15");
}

// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
// Cooling law (from metody_numeryczne_9 / _10):
// dT/dt = -alpha * T^4,   T(0) = 5903 K
// Exact: T(t) = T0*10000 / (1851232043943*t + 1e12)^(1/3)
// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
void test_cooling_law()
{
    cout << "\n[cooling law — physics ODE]:\n";
    const double alpha = 3e-12;
    const double T0 = 5903.0;
    const double t_end = 5903.0;

    auto f_cool = [alpha](double T, double) { return -alpha * T * T * T * T; };
    auto exact = [T0](double t) {
        return T0 * 10000.0 / pow(1851232043943.0 * t + 1e12, 1.0 / 3.0);
        };

    double T_rk4 = rk4_solve(f_cool, T0, t_end, 100000);
    double T_exact = exact(t_end);

    check(fabs(T_rk4 - T_exact) < 1.0, "cooling: RK4 error < 1 K at t=5903 s");

    double mseE = ode_mse(f_cool, exact, T0, t_end, 10000, 0);
    double mseRK = ode_mse(f_cool, exact, T0, t_end, 10000, 3);
    check(mseRK < mseE, "cooling: RK4 MSE < Euler MSE");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTS: ODE solvers\n";
    cout << "======================================\n";

    test_euler();
    test_heun();
    test_midpoint();
    test_rk4();
    test_accuracy_ordering();
    test_ode_mse();
    test_cooling_law();

    cout << "\n--------------------------------------\n";
    cout << "Results: " << passed << " passed, " << failed << " failed\n";
    return (failed == 0) ? 0 : 1;
}