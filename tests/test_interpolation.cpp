/**
 * @file test_interpolation.cpp
 * @brief Testy jednostkowe modułu interpolacji
 *
 * Kompilacja:  g++ -std=c++17 -I../include test_interpolation.cpp ../src/interpolation.cpp -o test_interp
 * Uruchomienie:    ./test_interp
 */

#include <iostream>
#include <cmath>
#include "interpolation.h"

using namespace std;

static int passed = 0, failed = 0;

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
// TEST 1: Lagrange odtwarza dokładny wielomian
// =========================================================================
void test_lagrange_poly()
{
    cout << "\n[lagrange_interpolate] dokladny na wezlach kwadratowych:\n";
    // f(x) = x^2; wezly w 0,1,2 -> wartosci 0,1,4
    vector<double> xs = { 0.0, 1.0, 2.0 };
    vector<double> ys = { 0.0, 1.0, 4.0 };

    // Musi być dokładny w węzłach
    check(near(lagrange_interpolate(0.0, xs, ys), 0.0), "dokladny w x=0");
    check(near(lagrange_interpolate(1.0, xs, ys), 1.0), "dokladny w x=1");
    check(near(lagrange_interpolate(2.0, xs, ys), 4.0), "dokladny w x=2");

    // Interpolacja w punkcie środkowym x^2 dla x=1.5 powinna dać dokładnie 2.25
    check(near(lagrange_interpolate(1.5, xs, ys), 2.25), "srodek x=1.5 -> 2.25");
}

// =========================================================================
// TEST 2: Lagrange – liniowy (2 węzły) jest dokładny w każdym punkcie
// =========================================================================
void test_lagrange_linear()
{
    cout << "\n[lagrange_interpolate] liniowy (2 wezly):\n";
    vector<double> xs = { 0.0, 4.0 };
    vector<double> ys = { 1.0, 5.0 };  // y = x + 1
    check(near(lagrange_interpolate(2.0, xs, ys), 3.0), "y=x+1 w x=2 -> 3");
    check(near(lagrange_interpolate(3.0, xs, ys), 4.0), "y=x+1 w x=3 -> 4");
}

// =========================================================================
// TEST 3: Różnice dzielone Newtona – ten sam wielomian kwadratowy co Lagrange
// =========================================================================
void test_newton_quadratic()
{
    cout << "\n[newton_interpolate] kwadratowy:\n";
    vector<double> xs = { 0.0, 1.0, 2.0 };
    vector<double> ys = { 0.0, 1.0, 4.0 };
    auto dd = newton_divided_differences(xs, ys);
    check(near(newton_interpolate(1.5, xs, dd), 2.25), "srodek x=1.5 -> 2.25");
    check(near(newton_interpolate(0.0, xs, dd), 0.0), "dokladny w x=0");
}

// =========================================================================
// TEST 4: Newton – sin(x) z 5 węzłami (dokładny w pobliżu środka)
// =========================================================================
void test_newton_sin()
{
    cout << "\n[newton_interpolate] sin(x) 5 wezlow:\n";
    vector<double> xs, ys;
    for (int i = 0; i <= 4; ++i) {
        double x = i * M_PI / 4.0;
        xs.push_back(x);
        ys.push_back(sin(x));
    }
    auto dd = newton_divided_differences(xs, ys);
    double val = newton_interpolate(M_PI / 6.0, xs, dd);  // sin(30°) = 0.5
    check(fabs(val - 0.5) < 0.01, "sin(pi/6) ~ 0.5 (tol 0.01)");
    check(near(newton_interpolate(0.0, xs, dd), 0.0, 1e-10), "dokladny w x=0");
}

// =========================================================================
// TEST 5: Schemat Hornera – znany wielomian 1 + 2x + 3x^2
// =========================================================================
void test_horner()
{
    cout << "\n[horner_eval]:\n";
    // p(x) = 1 + 2x + 3x^2
    vector<double> a = { 1.0, 2.0, 3.0 };
    check(near(horner_eval(a, 0.0), 1.0), "p(0) = 1");
    check(near(horner_eval(a, 1.0), 6.0), "p(1) = 6");
    check(near(horner_eval(a, -1.0), 2.0), "p(-1) = 2");
    check(near(horner_eval(a, 2.0), 17.0), "p(2) = 17");
}

// =========================================================================
// TEST 6: interpolation_mse – dokładne węzły => MSE == 0
// =========================================================================
void test_mse_exact()
{
    cout << "\n[interpolation_mse] zero przy interpolacji w wezlach:\n";
    vector<double> xs = { 0.0, 1.0, 2.0, 3.0 };
    vector<double> ys = { 0.0, 1.0, 4.0, 9.0 };  // x^2
    // Te same punkty jako dane i węzły
    double mse = interpolation_mse(xs, ys, xs, ys);
    check(mse < 1e-20, "MSE == 0 gdy wezly == punkty danych");

    // Tylko 2 węzły => wielomian nie dopasuje wszystkich 4 punktów => MSE > 0
    vector<double> xs2 = { xs[0], xs[3] };
    vector<double> ys2 = { ys[0], ys[3] };
    double mse2 = interpolation_mse(xs, ys, xs2, ys2);
    check(mse2 > 0.0, "MSE > 0 przy rzadszych wezlach");
}

// =========================================================================
int main()
{
    cout << "======================================\n";
    cout << "  TESTY: interpolacja\n";
    cout << "======================================\n";

    test_lagrange_poly();
    test_lagrange_linear();
    test_newton_quadratic();
    test_newton_sin();
    test_horner();
    test_mse_exact();

    cout << "\n--------------------------------------\n";
    cout << "Wyniki: " << passed << " zaliczone, " << failed << " niezaliczone\n";
    return (failed == 0) ? 0 : 1;
}