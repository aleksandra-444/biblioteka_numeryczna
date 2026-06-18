/**
 * @file test_integration.cpp
 * @brief Testy jednostkowe modułu całkowania numerycznego
 *
 * Kompilacja:  g++ -std=c++17 -I../include test_integration.cpp ../src/integration.cpp -o test_integ
 * Uruchomienie:    ./test_integ
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
// Całki referencyjne (dokładne):
//   int_0^1 x^2 dx = 1/3
//   int_0^pi sin(x) dx = 2
//   int_{-1}^{1} (1 - x^2) dx = 4/3
// =========================================================================

void test_trapezoidal()
{
    cout << "\n[trapezoidal]:\n";
    // int_0^1 x^2 dx = 1/3  (m=10000 -> bardzo dokladne)
    auto f1 = [](double x) { return x * x; };
    double I = trapezoidal(f1, 0.0, 1.0, 10000);
    check(fabs(I - 1.0 / 3.0) < 1e-6, "x^2 na [0,1] => 1/3  (tol 1e-6)");

    // int_0^pi sin(x) dx = 2
    auto f2 = [](double x) { return sin(x); };
    double I2 = trapezoidal(f2, 0.0, M_PI, 10000);
    check(fabs(I2 - 2.0) < 1e-6, "sin(x) na [0,pi] => 2  (tol 1e-6)");
}

void test_simpson()
{
    cout << "\n[simpson]:\n";
    // x^2: Simpson jest dokladny dla wielomianow do stopnia 3
    auto f1 = [](double x) { return x * x; };
    double I = simpson(f1, 0.0, 1.0, 100);
    check(fabs(I - 1.0 / 3.0) < 1e-12, "x^2 na [0,1] => 1/3  (tol 1e-12)");

    // sin(x)
    auto f2 = [](double x) { return sin(x); };
    double I2 = simpson(f2, 0.0, M_PI, 1000);
    check(fabs(I2 - 2.0) < 1e-10, "sin(x) na [0,pi] => 2  (tol 1e-10)");

    // Simpson jest dokladniejszy niz trapezy dla tego samego m
    double IT = trapezoidal(f2, 0.0, M_PI, 1000);
    double IS = simpson(f2, 0.0, M_PI, 1000);
    check(fabs(IS - 2.0) < fabs(IT - 2.0), "Simpson dokladniejszy niz Trapezy");
}

void test_gauss_legendre()
{
    cout << "\n[gauss_legendre / zlozony]:\n";
    // 1 - x^2 na [-1,1] = 4/3 (stopien 2 => dokladny dla n_nodes >= 2)
    auto f1 = [](double x) { return 1.0 - x * x; };
    double I = gauss_legendre(f1, -1.0, 1.0, 2);
    check(fabs(I - 4.0 / 3.0) < 1e-12, "1-x^2 na [-1,1] => 4/3, n=2, dokladny");

    // x^4 na [0,1] = 1/5; potrzeba n_nodes >= 3
    auto f2 = [](double x) { return pow(x, 4); };
    double I2 = gauss_legendre(f2, 0.0, 1.0, 3);
    check(fabs(I2 - 0.2) < 1e-10, "x^4 na [0,1] => 0.2, n=3");

    // Złożony: sin(x) na [0,pi], wiele podprzedziałów
    auto f3 = [](double x) { return sin(x); };
    double I3 = gauss_legendre_composite(f3, 0.0, M_PI, 100, 3);
    check(fabs(I3 - 2.0) < 1e-12, "sin(x) zlozony GL3 => 2");

    // Wiecej wezlow => wieksza dokladnosc dla gladkiej funkcji
    double I4a = gauss_legendre_composite(f3, 0.0, M_PI, 10, 2);
    double I4b = gauss_legendre_composite(f3, 0.0, M_PI, 10, 5);
    check(fabs(I4b - 2.0) < fabs(I4a - 2.0),
        "5 wezlow dokladniejsze niz 2 wezly (to samo m)");
}

void test_invalid_nodes()
{
    cout << "\n[gauss_legendre] nieprawidlowa liczba wezlow rzuca wyjatek:\n";
    bool threw = false;
    try {
        auto f = [](double x) { return x; };
        gauss_legendre(f, 0.0, 1.0, 7);  // 7 wezlow nie jest obsługiwane
    }
    catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "rzuca invalid_argument dla nieobslugiwanych n_nodes");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTY: calkowanie numeryczne\n";
    cout << "======================================\n";

    test_trapezoidal();
    test_simpson();
    test_gauss_legendre();
    test_invalid_nodes();

    cout << "\n--------------------------------------\n";
    cout << "Wyniki: " << passed << " zaliczone, " << failed << " niezaliczone\n";
    return (failed == 0) ? 0 : 1;
}