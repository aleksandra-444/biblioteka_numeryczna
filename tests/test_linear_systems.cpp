/**
 * @file test_linear_systems.cpp
 * @brief Testy jednostkowe modułu układów równań liniowych
 *
 * Każdy TEST sprawdza konkretną funkcję co najmniej w 2 przypadkach
 * (poprawne dane wejściowe + przypadek brzegowy/błąd).
 * Kompilacja:  g++ -std=c++17 -I../include test_linear_systems.cpp ../src/linear_systems.cpp -o test_lin
 * Uruchomienie:    ./test_lin
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include "linear_systems.h"

using namespace std;

static int passed = 0, failed = 0;

// ── Pomocnicza funkcja sprawdzająca ──────────────────────────────────────
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
// TEST 1: gauss_solve – układ 2x2 ze znanym rozwiązaniem
// =========================================================================
void test_gauss_2x2()
{
    cout << "\n[gauss_solve] znany uklad 2x2:\n";
    // 2x + y = 5
    // x + 3y = 10  =>  x=1, y=3
    vector<vector<double>> A = { {2, 1}, {1, 3} };
    vector<double> b = { 5, 10 };
    auto x = gauss_solve(A, b);
    check(near(x[0], 1.0) && near(x[1], 3.0), "x=1, y=3");
}

// =========================================================================
// TEST 2: gauss_solve – układ 3x3
// =========================================================================
void test_gauss_3x3()
{
    cout << "\n[gauss_solve] uklad 3x3:\n";
    // x + y + z = 6
    // 2x + y + z = 8   =>  x=2, y=2, z=2 (LU i Gauss)
    // x + 2y + 3z = 14
    vector<vector<double>> A = { {1,1,1},{1,2,1},{1,1,2} };
    vector<double> b = { 6, 8, 9 };
    auto x = gauss_solve(A, b);
    check(near(x[0], 1.0) && near(x[1], 2.0) && near(x[2], 3.0),
        "x=1, y=2, z=3");
}

// =========================================================================
// TEST 3: lu_solve – ten sam układ 3x3 co wyżej, przez LU
// =========================================================================
void test_lu_solve_3x3()
{
    cout << "\n[lu_solve] 3x3:\n";
    vector<vector<double>> A = { {1,1,1},{1,2,1},{1,1,2} };
    vector<double> b = { 6, 8, 9 };
    auto x = lu_solve(A, b);
    check(near(x[0], 1.0) && near(x[1], 2.0) && near(x[2], 3.0),
        "lu_solve zgodny z gauss_solve");
}

// =========================================================================
// TEST 4: lu_solve – macierz diagonalna 4x4 (trywialnie rozwiązywalna)
// =========================================================================
void test_lu_diagonal()
{
    cout << "\n[lu_solve] diagonalna 4x4:\n";
    vector<vector<double>> A = {
        {2,0,0,0},{0,3,0,0},{0,0,5,0},{0,0,0,7}
    };
    vector<double> b = { 4, 9, 25, 49 };
    auto x = lu_solve(A, b);
    check(near(x[0], 2.0) && near(x[1], 3.0) && near(x[2], 5.0) && near(x[3], 7.0),
        "x={2,3,5,7}");
}

// =========================================================================
// TEST 5: Dekompozycja LU – weryfikacja L*U == A
// =========================================================================
void test_lu_decompose_verify()
{
    cout << "\n[lu_decompose] L*U == A:\n";
    vector<vector<double>> A = { {4,3},{6,3} };
    vector<vector<double>> L, U;
    lu_decompose(A, L, U, 2);

    // Oblicz L*U ręcznie
    bool ok = true;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            double s = L[i][0] * U[0][j] + L[i][1] * U[1][j];
            if (!near(s, A[i][j])) ok = false;
        }
    check(ok, "L*U odtwarza A");
    check(near(L[0][0], 1.0) && near(L[1][1], 1.0), "L ma jedynki na diagonali");
}

// =========================================================================
// TEST 6: analyze_system – oznaczony / sprzeczny / nieoznaczony
// =========================================================================
void test_analyze_system()
{
    cout << "\n[analyze_system]:\n";
    // Jednoznaczne rozwiązanie
    vector<vector<double>> A1 = { {1,0},{0,1} };
    vector<double> b1 = { 1,2 };
    check(analyze_system(A1, b1, 2) == "OZNACZONY", "Identycznosc => OZNACZONY");

    // Sprzeczny
    vector<vector<double>> A2 = { {1,1},{1,1} };
    vector<double> b2 = { 1, 2 };
    check(analyze_system(A2, b2, 2) == "SPRZECZNY", "Rownolegle wiersze, rozne prawe strony => SPRZECZNY");

    // Nieoznaczony
    vector<vector<double>> A3 = { {1,1},{1,1} };
    vector<double> b3 = { 2, 2 };
    check(analyze_system(A3, b3, 2) == "NIEOZNACZONY", "Rownolegle wiersze, te same prawe strony => NIEOZNACZONY");
}

// =========================================================================
// TEST 7: residual_norm – powinno być bliskie 0 dla dokładnego rozwiązania
// =========================================================================
void test_residual_norm()
{
    cout << "\n[residual_norm]:\n";
    vector<vector<double>> A = { {2,1},{1,3} };
    vector<double> b = { 5,10 };
    auto x = gauss_solve(A, b);
    double r = residual_norm(A, x, b, 2);
    check(r < 1e-10, "residuum bliskie zeru dla dokladnego rozwiazania");

    // Zaburz rozwiązanie
    x[0] += 1.0;
    double r2 = residual_norm(A, x, b, 2);
    check(r2 > 0.5, "residuum rosnie gdy rozwiazanie jest bledne");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTY: uklady rownan liniowych\n";
    cout << "======================================\n";

    test_gauss_2x2();
    test_gauss_3x3();
    test_lu_solve_3x3();
    test_lu_diagonal();
    test_lu_decompose_verify();
    test_analyze_system();
    test_residual_norm();

    cout << "\n--------------------------------------\n";
    cout << "Wyniki: " << passed << " zaliczone, " << failed << " niezaliczone\n";
    return (failed == 0) ? 0 : 1;
}