/**
 * @file test_approximation.cpp
 * @brief Unit tests for approximation module
 *
 * Build:
 *   g++ -std=c++17 -I../include test_approximation.cpp \
 *       ../src/approximation.cpp ../src/linear_systems.cpp -o test_approx
 * Run:    ./test_approx
 */

#include <iostream>
#include <cmath>
#include "approximation.h"

using namespace std;

static int passed = 0, failed = 0;

static void check(bool cond, const string& name)
{
    if (cond) { cout << "  [PASS] " << name << "\n"; ++passed; }
    else { cout << "  [FAIL] " << name << "\n"; ++failed; }
}

// �������������������������������������������������������������������������
// TEST 1: approximating a constant is exact
// �������������������������������������������������������������������������
void test_constant()
{
    cout << "\n[least_squares_poly] approximate constant f(x)=5:\n";
    auto f = [](double) { return 5.0; };
    auto c = least_squares_poly(f, -1.0, 1.0, 0);  // degree 0
    check(fabs(eval_approx(c, 0.5) - 5.0) < 1e-8, "degree 0 constant => exact");
}

// �������������������������������������������������������������������������
// TEST 2: approximating a linear function with degree 1 is exact
// �������������������������������������������������������������������������
void test_linear_exact()
{
    cout << "\n[least_squares_poly] degree 1 approx of 2x+1:\n";
    auto f = [](double x) { return 2.0 * x + 1.0; };
    auto c = least_squares_poly(f, 0.0, 1.0, 1);
    check(fabs(eval_approx(c, 0.5) - 2.0) < 1e-6, "P(0.5) ~= 2.0");
    check(fabs(eval_approx(c, 0.0) - 1.0) < 1e-6, "P(0.0) ~= 1.0");
}

// �������������������������������������������������������������������������
// TEST 3: higher degree => lower MSE
// �������������������������������������������������������������������������
void test_degree_improves()
{
    cout << "\n[approx_mse] higher degree lowers error:\n";
    auto f = [](double x) { return (x - 1) * exp(sin(x * x)); };  // from lab 8
    auto c3 = least_squares_poly(f, -2.5, 1.0, 3);
    auto c10 = least_squares_poly(f, -2.5, 1.0, 10);
    double mse3 = approx_mse(f, c3, -2.5, 1.0);
    double mse10 = approx_mse(f, c10, -2.5, 1.0);
    check(mse10 < mse3, "degree 10 MSE < degree 3 MSE");

    double max3 = approx_max_error(f, c3, -2.5, 1.0);
    double max10 = approx_max_error(f, c10, -2.5, 1.0);
    check(max10 < max3, "degree 10 max error < degree 3 max error");
}

// �������������������������������������������������������������������������
// TEST 4: approximating x^2 with degree 2 must be nearly exact
// �������������������������������������������������������������������������
void test_quadratic_exact()
{
    cout << "\n[least_squares_poly] degree 2 approx of x^2:\n";
    auto f = [](double x) { return x * x; };
    auto c = least_squares_poly(f, -1.0, 1.0, 2);
    double mse = approx_mse(f, c, -1.0, 1.0);
    check(mse < 1e-20, "MSE ~ 0 for exact polynomial space");
}

// �������������������������������������������������������������������������
// TEST 5: eval_approx at degree-3 nodes (sanity check)
// �������������������������������������������������������������������������
void test_eval_approx()
{
    cout << "\n[eval_approx]:\n";
    // p(x) = 1 + 2x + 3x^2 + 4x^3, coeffs = {1,2,3,4}
    // via least_squares of this exact function
    auto f = [](double x) { return 1 + 2 * x + 3 * x * x + 4 * x * x * x; };
    auto c = least_squares_poly(f, -1.0, 1.0, 3);
    check(fabs(eval_approx(c, 0.0) - 1.0) < 1e-6, "P(0)=1");
    check(fabs(eval_approx(c, 1.0) - 10.0) < 1e-5, "P(1)=10");
    check(fabs(eval_approx(c, -1.0) - (-2.0)) < 1e-5, "P(-1)=-2");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTS: approximation\n";
    cout << "======================================\n";

    test_constant();
    test_linear_exact();
    test_degree_improves();
    test_quadratic_exact();
    test_eval_approx();

    cout << "\n--------------------------------------\n";
    cout << "Results: " << passed << " passed, " << failed << " failed\n";
    return (failed == 0) ? 0 : 1;
}