/**
 * @file example_roots_and_integrals.cpp
 * @brief Example 2: Root-finding + Numerical Integration
 *
 * Demonstrates all root-finding methods (bisection, regula falsi,
 * Newton, secant) and all integration methods (trapezoidal, Simpson,
 * Gauss-Legendre) on functions from laboratory exercises.
 *
 * Build:
 *   g++ -std=c++17 -I../include example_roots_and_integrals.cpp \
 *       ../src/nonlinear.cpp ../src/integration.cpp -o example_roots
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../include/nonlinear"
#include "../include/intergration"

using namespace std;

// ── Test functions ─────────────────────────────────────────────────────────
auto f1 = [](double x) { return x * x - 2.0; };        // root: sqrt(2)
auto df1 = [](double x) { return 2.0 * x; };

auto f2 = [](double x) { return cosh(x) - sqrt(x) - 1.0; };  // root ~1.2

// f(x) = (x+1)*exp(cos(x)), integral on [-14, 1] = -116.659583194914
auto f_int = [](double x) { return (x + 1.0) * exp(cos(x)); };
const double EXACT_INT = -116.659583194914;

// ── Print root-finding table ────────────────────────────────────────────────
static void print_root_table(const string& fname,
    double bisect_root, int bisect_it,
    double rf_root,
    double newt_root, int newt_it,
    double sec_root,
    double true_root)
{
    cout << "\n  Function: " << fname << "  (true root = " << true_root << ")\n";
    cout << "  " << string(70, '-') << "\n";
    cout << setw(20) << "Method"
        << setw(18) << "Root"
        << setw(16) << "|error|"
        << setw(12) << "Iters" << "\n";
    cout << "  " << string(70, '-') << "\n";
    auto row = [&](const string& m, double r, int it) {
        cout << "  " << setw(18) << m
            << setw(18) << fixed << setprecision(12) << r
            << setw(16) << fabs(r - true_root)
            << setw(12) << it << "\n";
        };
    row("Bisection", bisect_root, bisect_it);
    row("Regula Falsi", rf_root, -1);
    row("Newton", newt_root, newt_it);
    row("Secant", sec_root, -1);
    cout << "  " << string(70, '-') << "\n";
}

int main()
{
    cout << "==========================================================\n";
    cout << "  EXAMPLE 2: Root-Finding & Numerical Integration\n";
    cout << "==========================================================\n";

    cout << fixed << setprecision(12);

    // ════════════════════════════════════════════════════════════════
    // PART 1: Root finding
    // ════════════════════════════════════════════════════════════════
    cout << "\n======= ROOT FINDING =======\n";

    double b1 = bisection(f1, 1.0, 2.0);
    double rf1 = regula_falsi(f1, 1.0, 2.0);
    double n1 = newton_method(f1, df1, 2.0);
    double s1 = secant_method(f1, 2.0, 1.5);
    int ib1 = bisection_iterations(f1, 1.0, 2.0);
    int in1 = newton_iterations(f1, df1, 2.0);

    print_root_table("x^2 - 2 = 0", b1, ib1, rf1, n1, in1, s1, sqrt(2.0));

    double b2 = bisection(f2, 0.5, 1.5);
    double rf2 = regula_falsi(f2, 0.5, 1.5);
    // Numerical Newton (no analytical derivative available)
    double n2 = newton_numeric(f2, 1.0);
    double s2 = secant_method(f2, 0.5, 1.5);
    int ib2 = bisection_iterations(f2, 0.5, 1.5);
    // true root of cosh(x)-sqrt(x)-1 near 1.2 (no closed form)
    cout << "\n  Function: cosh(x) - sqrt(x) - 1 = 0\n";
    cout << "  " << string(60, '-') << "\n";
    cout << setw(20) << "Bisection:" << setw(20) << b2 << "  iters=" << ib2 << "\n";
    cout << setw(20) << "Regula Falsi:" << setw(20) << rf2 << "\n";
    cout << setw(20) << "Newton(num):" << setw(20) << n2 << "\n";
    cout << setw(20) << "Secant:" << setw(20) << s2 << "\n";
    cout << "  All methods agree: |max-min| = "
        << fabs(max({ b2,rf2,n2,s2 }) - min({ b2,rf2,n2,s2 })) << "\n";

    // ════════════════════════════════════════════════════════════════
    // PART 2: Numerical integration of (x+1)*exp(cos(x)) on [-14,1]
    // ════════════════════════════════════════════════════════════════
    cout << "\n\n======= NUMERICAL INTEGRATION =======\n";
    cout << "  f(x) = (x+1)*exp(cos(x)),  interval [-14, 1]\n";
    cout << "  Exact value = " << EXACT_INT << "\n\n";

    auto print_int = [&](const string& name, double val) {
        cout << "  " << setw(30) << name << " = " << setw(18) << val
            << "   |error| = " << fabs(val - EXACT_INT) << "\n";
        };

    for (int m : {100, 1000, 10000}) {
        cout << "  --- m = " << m << " sub-intervals ---\n";
        print_int("Trapezoidal", trapezoidal(f_int, -14.0, 1.0, m));
        print_int("Simpson", simpson(f_int, -14.0, 1.0, m));
        print_int("Gauss-Legendre n=2", gauss_legendre_composite(f_int, -14.0, 1.0, m, 2));
        print_int("Gauss-Legendre n=3", gauss_legendre_composite(f_int, -14.0, 1.0, m, 3));
        print_int("Gauss-Legendre n=4", gauss_legendre_composite(f_int, -14.0, 1.0, m, 4));
        cout << "\n";
    }

    // ── Save convergence to CSV ──────────────────────────────────────
    ofstream csv("integration_convergence.csv");
    csv << "m;Trapezoidal;Simpson;Gauss2;Gauss3;Gauss4\n";
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
    cout << "  Convergence data saved to: integration_convergence.csv\n";

    return 0;
}