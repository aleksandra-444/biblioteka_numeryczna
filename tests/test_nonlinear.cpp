/**
 * @file test_nonlinear.cpp
 * @brief Testy jednostkowe modułu znajdowania pierwiastków równań nieliniowych
 *
 * Kompilacja:  g++ -std=c++17 -I../include test_nonlinear.cpp ../src/nonlinear.cpp -o test_nl
 * Uruchomienie:    ./test_nl
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

// ─────────────────────────────────────────────────────────────────────────
// Funkcje testowe
// ─────────────────────────────────────────────────────────────────────────
auto f1 = [](double x) { return x * x - 2.0; };        // pierwiastek: sqrt(2) ~ 1.41421
auto df1 = [](double x) { return 2.0 * x; };

auto f2 = [](double x) { return x * x * x + x * x - 3 * x - 3; }; // pierwiastki: ±sqrt(3), -1
auto df2 = [](double x) { return 3 * x * x + 2 * x - 3; };

auto f3 = [](double x) { return cosh(x) - sqrt(x) - 1.0; }; // pierwiastek w pobliżu 1.2
auto f4 = [](double x) { return pow(x, 5) - 5 * pow(x, 4) + 10 * pow(x, 3)
- 10 * x * x + 5 * x - 1.0; };  // pierwiastek w x=1 (krotność 5)

// =========================================================================
void test_bisection()
{
    cout << "\n[bisection]:\n";
    double r = bisection(f1, 1.0, 2.0);
    check(fabs(r - sqrt(2.0)) < 1e-9, "sqrt(2): blad bisekcji < 1e-9");

    double r2 = bisection(f2, 1.0, 2.0);
    check(fabs(r2 - sqrt(3.0)) < 1e-9, "sqrt(3): blad bisekcji < 1e-9");

    // Błędny przedział => wyjątek
    bool threw = false;
    try { bisection(f1, 2.0, 3.0); }
    catch (const invalid_argument&) { threw = true; }
    check(threw, "throws when f(a)*f(b) >= 0");

    // Liczba iteracji powinna być ograniczona
    int iters = bisection_iterations(f1, 1.0, 2.0);
    check(iters > 0 && iters <= NL_MAX_ITER, "liczba iteracji w [1, MAX]");
}

// =========================================================================
void test_regula_falsi()
{
    cout << "\n[regula_falsi]:\n";
    double r = regula_falsi(f1, 1.0, 2.0);
    check(fabs(r - sqrt(2.0)) < 1e-9, "sqrt(2): blad regula falsi < 1e-9");

    double r2 = regula_falsi(f3, 0.5, 1.5);
    check(fabs(f3(r2)) < 1e-8, "cosh-sqrt-1: |f(pierwiastek)| < 1e-8");

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
    check(fabs(r - sqrt(2.0)) < 1e-12, "sqrt(2): blad Newtona < 1e-12");

    double r_neg = newton_method(f1, df1, -2.0);
    check(fabs(r_neg + sqrt(2.0)) < 1e-12, "pierwiastek ujemny: -sqrt(2), blad < 1e-12");

    // Pierwiastki wielomianu
    double r2 = newton_method(f2, df2, 2.0);
    check(fabs(r2 - sqrt(3.0)) < 1e-10, "Wielomian sqrt(3), blad < 1e-10");

    int it = newton_iterations(f1, df1, 2.0);
    int it_bis = bisection_iterations(f1, 1.0, 2.0);
    check(it < it_bis, "Newton zbiega szybciej niz Bisekcja");
}

// =========================================================================
void test_newton_numeric()
{
    cout << "\n[newton_numeric]:\n";
    double r = newton_numeric(f1, 2.0);
    check(fabs(r - sqrt(2.0)) < 1e-10, "sqrt(2): blad Newtona numerycznego < 1e-10");

    // Powinien dawać ten sam wynik co Newton analityczny
    double r_anal = newton_method(f1, df1, 2.0);
    check(fabs(r - r_anal) < 1e-9, "Newton numeryczny ~ Newton analityczny");
}

// =========================================================================
void test_secant()
{
    cout << "\n[secant_method]:\n";
    double r = secant_method(f1, 2.0, 1.5);
    check(fabs(r - sqrt(2.0)) < 1e-10, "sqrt(2): blad metody siecznych < 1e-10");

    double r2 = secant_method(f2, -0.5, -1.2);
    check(fabs(r2 + 1.0) < 1e-9, "Pierwiastek wielomianu -1, blad < 1e-9");
}

// =========================================================================
// Przypadek szczególny: wielomian z pierwiastkiem 5-krotnym w x=1
// Wszystkie metody powinny zbiegać, choć wolniej
// =========================================================================
void test_multiple_root()
{
    cout << "\n[pierwiastek wielokrotny x=1 (krotnosc 5)]:\n";
    auto df4 = [](double x) { return 5 * pow(x, 4) - 20 * pow(x, 3) + 30 * x * x - 20 * x + 5; };
    double r = bisection(f4, 0.9, 1.1, 1e-9);
    check(fabs(r - 1.0) < 1e-6, "bisekcja: pierwiastek bliski 1.0, tol 1e-6");

    double r2 = newton_method(f4, df4, 0.9, 1e-9);
    check(fabs(r2 - 1.0) < 1e-3, "newton: pierwiastek bliski 1.0 (wielokrotny, tol 1e-3)");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTY: rownienia nieliniowe\n";
    cout << "======================================\n";

    test_bisection();
    test_regula_falsi();
    test_newton();
    test_newton_numeric();
    test_secant();
    test_multiple_root();

    cout << "\n--------------------------------------\n";
    cout << "Wyniki: " << passed << " zaliczone, " << failed << " niezaliczone\n";
    return (failed == 0) ? 0 : 1;
}