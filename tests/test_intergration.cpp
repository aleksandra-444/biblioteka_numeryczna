/**
 * @file test_integration.cpp
 * @brief Unit tests for integration module
 *
 * Build:  g++ -std=c++17 -I../include test_integration.cpp ../src/integration.cpp -o test_integ
 * Run:    ./test_integ
 */

#include <iostream>
#include <cmath>
#include "integration.h"

using namespace std;

static int passed = 0, failed = 0;

static void check(bool cond, const string& name)
{
    if (cond) { cout << "  [PASS] " << name << "\n"; ++passed; }
    else { cout << "  [FAIL] " << name << "\n"; ++failed; }
}

// =========================================================================
// Reference integrals (exact):
//   int_0^1 x^2 dx = 1/3
//   int_0^pi sin(x) dx = 2
//   int_{-1}^{1} (1 - x^2) dx = 4/3
// =========================================================================

void test_trapezoidal()
{
    cout << "\n[trapezoidal]:\n";
    // int_0^1 x^2 dx = 1/3  (m=10000 -> very accurate)
    auto f1 = [](double x) { return x * x; };
    double I = trapezoidal(f1, 0.0, 1.0, 10000);
    check(fabs(I - 1.0 / 3.0) < 1e-6, "x^2 on [0,1] => 1/3  (tol 1e-6)");

    // int_0^pi sin(x) dx = 2
    auto f2 = [](double x) { return sin(x); };
    double I2 = trapezoidal(f2, 0.0, M_PI, 10000);
    check(fabs(I2 - 2.0) < 1e-6, "sin(x) on [0,pi] => 2  (tol 1e-6)");
}

void test_simpson()
{
    cout << "\n[simpson]:\n";
    // x^2: Simpson is exact for polynomials up to degree 3
    auto f1 = [](double x) { return x * x; };
    double I = simpson(f1, 0.0, 1.0, 100);
    check(fabs(I - 1.0 / 3.0) < 1e-12, "x^2 on [0,1] => 1/3  (tol 1e-12)");

    // sin(x)
    auto f2 = [](double x) { return sin(x); };
    double I2 = simpson(f2, 0.0, M_PI, 1000);
    check(fabs(I2 - 2.0) < 1e-10, "sin(x) on [0,pi] => 2  (tol 1e-10)");

    // Simpson is more accurate than trapezoidal for same m
    double IT = trapezoidal(f2, 0.0, M_PI, 1000);
    double IS = simpson(f2, 0.0, M_PI, 1000);
    check(fabs(IS - 2.0) < fabs(IT - 2.0), "Simpson more accurate than Trapz");
}

void test_gauss_legendre()
{
    cout << "\n[gauss_legendre / composite]:\n";
    // 1 - x^2 on [-1,1] = 4/3 (degree 2 => exact with n_nodes >= 2)
    auto f1 = [](double x) { return 1.0 - x * x; };
    double I = gauss_legendre(f1, -1.0, 1.0, 2);
    check(fabs(I - 4.0 / 3.0) < 1e-12, "1-x^2 on [-1,1] => 4/3, n=2, exact");

    // x^4 on [0,1] = 1/5; need n_nodes >= 3
    auto f2 = [](double x) { return pow(x, 4); };
    double I2 = gauss_legendre(f2, 0.0, 1.0, 3);
    check(fabs(I2 - 0.2) < 1e-10, "x^4 on [0,1] => 0.2, n=3");

    // Composite: sin(x) on [0,pi], many sub-intervals
    auto f3 = [](double x) { return sin(x); };
    double I3 = gauss_legendre_composite(f3, 0.0, M_PI, 100, 3);
    check(fabs(I3 - 2.0) < 1e-12, "sin(x) composite GL3 => 2");

    // More nodes => higher accuracy for smooth function
    double I4a = gauss_legendre_composite(f3, 0.0, M_PI, 10, 2);
    double I4b = gauss_legendre_composite(f3, 0.0, M_PI, 10, 5);
    check(fabs(I4b - 2.0) < fabs(I4a - 2.0),
        "5 nodes more accurate than 2 nodes (same m)");
}

void test_invalid_nodes()
{
    cout << "\n[gauss_legendre] invalid n_nodes throws:\n";
    bool threw = false;
    try {
        auto f = [](double x) { return x; };
        gauss_legendre(f, 0.0, 1.0, 7);  // 7 nodes not supported
    }
    catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "throws invalid_argument for unsupported n_nodes");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTS: integration\n";
    cout << "======================================\n";

    test_trapezoidal();
    test_simpson();
    test_gauss_legendre();
    test_invalid_nodes();

    cout << "\n--------------------------------------\n";
    cout << "Results: " << passed << " passed, " << failed << " failed\n";
    return (failed == 0) ? 0 : 1;
}