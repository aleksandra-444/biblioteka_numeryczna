/**
 * @file example_roots_and_integrals.cpp
 * @brief Przykład 2: Znajdowanie pierwiastków + Całkowanie numeryczne
 *
 * Demonstracja wszystkich metod znajdowania pierwiastków (bisekcja, regula falsi,
 * Newton, sieczna) oraz wszystkich metod całkowania (trapezy, Simpson,
 * Gauss-Legendre) na funkcjach z ćwiczeń laboratoryjnych.
 *
 * Kompilacja:
 *   g++ -std=c++17 -I../include example_roots_and_integrals.cpp \
 *       ../src/nonlinear.cpp ../src/integration.cpp -o example_roots
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../include/nonlinear.h"
#include "../include/integration.h"

using namespace std;

// ── Funkcje testowe ───────────────────────────────────────────────────────
auto f1 = [](double x) { return x * x - 2.0; };        // pierwiastek: sqrt(2)
auto df1 = [](double x) { return 2.0 * x; };

auto f2 = [](double x) { return cosh(x) - sqrt(x) - 1.0; };  // pierwiastek ~1.2

// f(x) = (x+1)*exp(cos(x)), całka na [-14, 1] = -116.659583194914
auto f_int = [](double x) { return (x + 1.0) * exp(cos(x)); };
const double EXACT_INT = -116.659583194914;

// ── Wypisz tabelę metod znajdowania pierwiastków ──────────────────────────
static void print_root_table(const string& fname,
    double bisect_root, int bisect_it,
    double rf_root,
    double newt_root, int newt_it,
    double sec_root,
    double true_root)
{
    cout << "\n  Funkcja: " << fname << "  (dokladny pierwiastek = " << true_root << ")\n";
    cout << "  " << string(70, '-') << "\n";
    cout << setw(20) << "Metoda"
        << setw(18) << "Pierwiastek"
        << setw(16) << "|blad|"
        << setw(12) << "Iter." << "\n";
    cout << "  " << string(70, '-') << "\n";
    auto row = [&](const string& m, double r, int it) {
        cout << "  " << setw(18) << m
            << setw(18) << fixed << setprecision(12) << r
            << setw(16) << fabs(r - true_root)
            << setw(12) << it << "\n";
        };
    row("Bisekcja", bisect_root, bisect_it);
    row("Regula Falsi", rf_root, -1);
    row("Newton", newt_root, newt_it);
    row("Sieczna", sec_root, -1);
    cout << "  " << string(70, '-') << "\n";
}

int main()
{
    cout << "==========================================================\n";
    cout << "  PRZYKLAD 2: Znajdowanie pierwiastkow i calkowanie\n";
    cout << "==========================================================\n";

    cout << fixed << setprecision(12);

    // ════════════════════════════════════════════════════════════════
    // CZESC 1: Znajdowanie pierwiastków
    // ════════════════════════════════════════════════════════════════
    cout << "\n======= ZNAJDOWANIE PIERWIASTKOW =======\n";

    double b1 = bisection(f1, 1.0, 2.0);
    double rf1 = regula_falsi(f1, 1.0, 2.0);
    double n1 = newton_method(f1, df1, 2.0);
    double s1 = secant_method(f1, 2.0, 1.5);
    int ib1 = bisection_iterations(f1, 1.0, 2.0);
    int in1 = newton_iterations(f1, df1, 2.0);

    print_root_table("x^2 - 2 = 0", b1, ib1, rf1, n1, in1, s1, sqrt(2.0));

    double b2 = bisection(f2, 0.5, 1.5);
    double rf2 = regula_falsi(f2, 0.5, 1.5);
    // Newton numeryczny (brak analitycznej pochodnej)
    double n2 = newton_numeric(f2, 1.0);
    double s2 = secant_method(f2, 0.5, 1.5);
    int ib2 = bisection_iterations(f2, 0.5, 1.5);
    // dokladny pierwiastek cosh(x)-sqrt(x)-1 w pobliżu 1.2 (brak postaci zamknietej)
    cout << "\n  Funkcja: cosh(x) - sqrt(x) - 1 = 0\n";
    cout << "  " << string(60, '-') << "\n";
    cout << setw(20) << "Bisekcja:" << setw(20) << b2 << "  iter=" << ib2 << "\n";
    cout << setw(20) << "Regula Falsi:" << setw(20) << rf2 << "\n";
    cout << setw(20) << "Newton(num):" << setw(20) << n2 << "\n";
    cout << setw(20) << "Sieczna:" << setw(20) << s2 << "\n";
    cout << "  Wszystkie metody zgodne: |max-min| = "
        << fabs(max({ b2,rf2,n2,s2 }) - min({ b2,rf2,n2,s2 })) << "\n";

    // ════════════════════════════════════════════════════════════════
    // CZESC 2: Całkowanie numeryczne (x+1)*exp(cos(x)) na [-14,1]
    // ════════════════════════════════════════════════════════════════
    cout << "\n\n======= CALKOWANIE NUMERYCZNE =======\n";
    cout << "  f(x) = (x+1)*exp(cos(x)),  przedzial [-14, 1]\n";
    cout << "  Wartosc dokladna = " << EXACT_INT << "\n\n";

    auto print_int = [&](const string& name, double val) {
        cout << "  " << setw(30) << name << " = " << setw(18) << val
            << "   |blad| = " << fabs(val - EXACT_INT) << "\n";
        };

    for (int m : {100, 1000, 10000}) {
        cout << "  --- m = " << m << " podprzedzialy ---\n";
        print_int("Trapezy", trapezoidal(f_int, -14.0, 1.0, m));
        print_int("Simpson", simpson(f_int, -14.0, 1.0, m));
        print_int("Gauss-Legendre n=2", gauss_legendre_composite(f_int, -14.0, 1.0, m, 2));
        print_int("Gauss-Legendre n=3", gauss_legendre_composite(f_int, -14.0, 1.0, m, 3));
        print_int("Gauss-Legendre n=4", gauss_legendre_composite(f_int, -14.0, 1.0, m, 4));
        cout << "\n";
    }

    // ── Zapis zbieżności do CSV ───────────────────────────────────────
    ofstream csv("integration_convergence.csv");
    csv << "m;Trapezy;Simpson;Gauss2;Gauss3;Gauss4\n";
    for (int m = 10; m <= 10000; m = (m < 100 ? m + 10 : m < 1000 ? m + 100 : m + 1000)) {
        double IT = trapezoidal(f_int, -14.0, 1.0, m);
        double IS = simpson(f_int, -14.0, 1.0, m);
        double G2 = gauss_legendre_composite(f_int, -14.0, 1.0, m, 2);
        double G3 = gauss_legendre_composite(f_int, -14.0, 1.0, m, 3);
        double G4 = gauss_legendre_composite(f_int, -14.0, 1.0, m, 4);
        csv << m << ";" << fabs(IT - EXACT_INT) << ";" << fabs(IS - EXACT_INT)
            << ";" << fabs(G2 - EXACT_INT) << ";" << fabs(G3 - EXACT_INT)
            << ";" << fabs(G4 - EXACT_INT) << "\n";
    }
    csv.close();
    cout << "  Dane zbieznosci zapisano do: integration_convergence.csv\n";

    return 0;
}