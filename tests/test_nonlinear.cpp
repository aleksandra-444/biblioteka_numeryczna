/**
 * @file test_nonlinear.cpp
 * @brief Unit tests for nonlinear root-finding module
 *
 * Build:  g++ -std=c++17 -I../include test_nonlinear.cpp ../src/nonlinear.cpp -o test_nl
 * Run:    ./test_nl
 */

#include <iostream>
#include <cmath>
#include <stdexcept>
#include "nonlinear.h"

using namespace std;

static int passed = 0, failed = 0;

static void check(bool cond, const string& name)
{
    if (cond) { cout << "  [PASS] " << name << "\n"; ++passed; }
    else { cout << "  [FAIL] " << name << "\n"; ++failed; }
}

// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
// Test functions (same as in metody_numeryczne_11 / _12)
// ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
auto f1 = [](double x) { return x * x - 2.0; };        // root: sqrt(2) ~ 1.41421
auto df1 = [](double x) { return 2.0 * x; };

auto f2 = [](double x) { return x * x * x + x * x - 3 * x - 3; }; // roots: ±sqrt(3), -1
auto df2 = [](double x) { return 3 * x * x + 2 * x - 3; };

auto f3 = [](double x) { return cosh(x) - sqrt(x) - 1.0; }; // root near 1.2
auto f4 = [](double x) { return pow(x, 5) - 5 * pow(x, 4) + 10 * pow(x, 3)
- 10 * x * x + 5 * x - 1.0; };  // root at x=1 (multiplicity 5)

// =========================================================================
void test_bisection()
{
    cout << "\n[bisection]:\n";
    double r = bisection(f1, 1.0, 2.0);
    check(fabs(r - sqrt(2.0)) < 1e-9, "sqrt(2): bisection error < 1e-9");

    double r2 = bisection(f2, 1.0, 2.0);
    check(fabs(r2 - sqrt(3.0)) < 1e-9, "sqrt(3): bisection error < 1e-9");

    // Wrong bracket => throws
    bool threw = false;
    try { bisection(f1, 2.0, 3.0); }
    catch (const invalid_argument&) { threw = true; }
    check(threw, "throws when f(a)*f(b) >= 0");

    // Iteration count should be bounded
    int iters = bisection_iterations(f1, 1.0, 2.0);
    check(iters > 0 && iters <= NL_MAX_ITER, "iteration count in [1, MAX]");
}

// =========================================================================
void test_regula_falsi()
{
    cout << "\n[regula_falsi]:\n";
    double r = regula_falsi(f1, 1.0, 2.0);
    check(fabs(r - sqrt(2.0)) < 1e-9, "sqrt(2): regula falsi error < 1e-9");

    double r2 = regula_falsi(f3, 0.5, 1.5);
    check(fabs(f3(r2)) < 1e-8, "cosh-sqrt-1: |f(root)| < 1e-8");

    bool threw = false;
    try { regula_falsi(f1, 2.0, 3.0); }
    catch (const invalid_argument&) { threw = true; }
    check(threw, "throws when f(a)*f(b) >= 0");
}

// =========================================================================
void test_newton()
{
    cout << "\n[newton_method]:\n";
    double r = newton_method(f1, df1, 2.0);
    check(fabs(r - sqrt(2.0)) < 1e-12, "sqrt(2): Newton error < 1e-12");

    double r_neg = newton_method(f1, df1, -2.0);
    check(fabs(r_neg + sqrt(2.0)) < 1e-12, "negative root: -sqrt(2), error < 1e-12");

    // Polynomial roots
    double r2 = newton_method(f2, df2, 2.0);
    check(fabs(r2 - sqrt(3.0)) < 1e-10, "Polynomial sqrt(3), error < 1e-10");

    int it = newton_iterations(f1, df1, 2.0);
    int it_bis = bisection_iterations(f1, 1.0, 2.0);
    check(it < it_bis, "Newton converges faster than Bisection");
}

// =========================================================================
void test_newton_numeric()
{
    cout << "\n[newton_numeric]:\n";
    double r = newton_numeric(f1, 2.0);
    check(fabs(r - sqrt(2.0)) < 1e-10, "sqrt(2): numeric Newton error < 1e-10");

    // Should give same answer as analytical Newton
    double r_anal = newton_method(f1, df1, 2.0);
    check(fabs(r - r_anal) < 1e-9, "numeric ~ analytic Newton result");
}

// =========================================================================
void test_secant()
{
    cout << "\n[secant_method]:\n";
    double r = secant_method(f1, 2.0, 1.5);
    check(fabs(r - sqrt(2.0)) < 1e-10, "sqrt(2): secant error < 1e-10");

    double r2 = secant_method(f2, -0.5, -1.2);
    check(fabs(r2 + 1.0) < 1e-9, "Polynomial root -1, error < 1e-9");
}

// =========================================================================
// Degenerate case: polynomial with 5-fold root at x=1
// All methods should converge but slowly
// =========================================================================
void test_multiple_root()
{
    cout << "\n[multiple root x=1 (multiplicity 5)]:\n";
    auto df4 = [](double x) { return 5 * pow(x, 4) - 20 * pow(x, 3) + 30 * x * x - 20 * x + 5; };
    double r = bisection(f4, 0.9, 1.1, 1e-9);
    check(fabs(r - 1.0) < 1e-6, "bisection: root near 1.0, tol 1e-6");

    double r2 = newton_method(f4, df4, 0.9, 1e-9);
    check(fabs(r2 - 1.0) < 1e-3, "newton: root near 1.0 (multiple root, tol 1e-3)");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTS: nonlinear root-finding\n";
    cout << "======================================\n";

    test_bisection();
    test_regula_falsi();
    test_newton();
    test_newton_numeric();
    test_secant();
    test_multiple_root();

    cout << "\n--------------------------------------\n";
    cout << "Results: " << passed << " passed, " << failed << " failed\n";
    return (failed == 0) ? 0 : 1;
}