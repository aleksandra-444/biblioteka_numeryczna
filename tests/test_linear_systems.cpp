/**
 * @file test_linear_systems.cpp
 * @brief Unit tests for linear_systems module
 *
 * Each TEST checks a specific function with at least 2 cases
 * (correct input + edge/error case).
 * Build:  g++ -std=c++17 -I../include test_linear_systems.cpp ../src/linear_systems.cpp -o test_lin
 * Run:    ./test_lin
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include "linear_systems.h"

using namespace std;

static int passed = 0, failed = 0;

// ¦¦ Tiny assert helper ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
static void check(bool cond, const string& name)
{
    if (cond) { cout << "  [PASS] " << name << "\n"; ++passed; }
    else { cout << "  [FAIL] " << name << "\n"; ++failed; }
}

static bool near(double a, double b, double tol = 1e-7)
{
    return fabs(a - b) < tol;
}

// =========================================================================
// TEST 1: gauss_solve — 2x2 system with known solution
// =========================================================================
void test_gauss_2x2()
{
    cout << "\n[gauss_solve] 2x2 known system:\n";
    // 2x + y = 5
    // x + 3y = 10  =>  x=1, y=3
    vector<vector<double>> A = { {2, 1}, {1, 3} };
    vector<double> b = { 5, 10 };
    auto x = gauss_solve(A, b);
    check(near(x[0], 1.0) && near(x[1], 3.0), "x=1, y=3");
}

// =========================================================================
// TEST 2: gauss_solve — 3x3 system
// =========================================================================
void test_gauss_3x3()
{
    cout << "\n[gauss_solve] 3x3 system:\n";
    // x + y + z = 6
    // 2x + y + z = 8   =>  x=2, y=2, z=2 (LU and Gauss)
    // x + 2y + 3z = 14
    vector<vector<double>> A = { {1,1,1},{1,2,1},{1,1,2} };
    vector<double> b = { 6, 8, 9 };
    auto x = gauss_solve(A, b);
    check(near(x[0], 1.0) && near(x[1], 2.0) && near(x[2], 3.0),
        "x=1, y=2, z=3");
}

// =========================================================================
// TEST 3: lu_solve — same 3x3 as above, using LU
// =========================================================================
void test_lu_solve_3x3()
{
    cout << "\n[lu_solve] 3x3:\n";
    vector<vector<double>> A = { {1,1,1},{1,2,1},{1,1,2} };
    vector<double> b = { 6, 8, 9 };
    auto x = lu_solve(A, b);
    check(near(x[0], 1.0) && near(x[1], 2.0) && near(x[2], 3.0),
        "lu_solve matches gauss_solve result");
}

// =========================================================================
// TEST 4: lu_solve — 4x4 diagonal matrix (trivially solvable)
// =========================================================================
void test_lu_diagonal()
{
    cout << "\n[lu_solve] 4x4 diagonal:\n";
    vector<vector<double>> A = {
        {2,0,0,0},{0,3,0,0},{0,0,5,0},{0,0,0,7}
    };
    vector<double> b = { 4, 9, 25, 49 };
    auto x = lu_solve(A, b);
    check(near(x[0], 2.0) && near(x[1], 3.0) && near(x[2], 5.0) && near(x[3], 7.0),
        "x={2,3,5,7}");
}

// =========================================================================
// TEST 5: LU decomposition — verify L*U == A
// =========================================================================
void test_lu_decompose_verify()
{
    cout << "\n[lu_decompose] L*U == A:\n";
    vector<vector<double>> A = { {4,3},{6,3} };
    vector<vector<double>> L, U;
    lu_decompose(A, L, U, 2);

    // Compute L*U manually
    bool ok = true;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            double s = L[i][0] * U[0][j] + L[i][1] * U[1][j];
            if (!near(s, A[i][j])) ok = false;
        }
    check(ok, "L*U reconstructs A");
    check(near(L[0][0], 1.0) && near(L[1][1], 1.0), "L has unit diagonal");
}

// =========================================================================
// TEST 6: analyze_system — definite / inconsistent / underdetermined
// =========================================================================
void test_analyze_system()
{
    cout << "\n[analyze_system]:\n";
    // Unique solution
    vector<vector<double>> A1 = { {1,0},{0,1} };
    vector<double> b1 = { 1,2 };
    check(analyze_system(A1, b1, 2) == "OZNACZONY", "Identity => OZNACZONY");

    // Inconsistent
    vector<vector<double>> A2 = { {1,1},{1,1} };
    vector<double> b2 = { 1, 2 };
    check(analyze_system(A2, b2, 2) == "SPRZECZNY", "Parallel rows, diff RHS => SPRZECZNY");

    // Underdetermined
    vector<vector<double>> A3 = { {1,1},{1,1} };
    vector<double> b3 = { 2, 2 };
    check(analyze_system(A3, b3, 2) == "NIEOZNACZONY", "Parallel rows, same RHS => NIEOZNACZONY");
}

// =========================================================================
// TEST 7: residual_norm — should be near 0 for exact solution
// =========================================================================
void test_residual_norm()
{
    cout << "\n[residual_norm]:\n";
    vector<vector<double>> A = { {2,1},{1,3} };
    vector<double> b = { 5,10 };
    auto x = gauss_solve(A, b);
    double r = residual_norm(A, x, b, 2);
    check(r < 1e-10, "residual near zero for exact solution");

    // Perturb the solution
    x[0] += 1.0;
    double r2 = residual_norm(A, x, b, 2);
    check(r2 > 0.5, "residual grows when solution is wrong");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTS: linear_systems\n";
    cout << "======================================\n";

    test_gauss_2x2();
    test_gauss_3x3();
    test_lu_solve_3x3();
    test_lu_diagonal();
    test_lu_decompose_verify();
    test_analyze_system();
    test_residual_norm();

    cout << "\n--------------------------------------\n";
    cout << "Results: " << passed << " passed, " << failed << " failed\n";
    return (failed == 0) ? 0 : 1;
}