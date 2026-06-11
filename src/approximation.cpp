#include "../include/approximation.h"
#include "../include/linear_systems.h"   // for gauss_solve
#include <cmath>
using std::vector;
using std::function;

// basis: phi_i(x) = x^i
static double phi(int i, double x) { return std::pow(x, i); }

// Simpson integration helper (always uses even number of steps)
static double simps_int(function<double(double)> g, double a, double b, int n)
{
    if (n % 2 != 0) ++n;
    double h = (b - a) / n;
    double s = g(a) + g(b);
    for (int k = 1; k < n; ++k)
        s += (k % 2 == 0 ? 2.0 : 4.0) * g(a + k * h);
    return (h / 3.0) * s;
}

// ── Build and solve normal equations ─────────────────────────────────────
vector<double> least_squares_poly(function<double(double)> f,
    double a, double b, int degree, int n_int)
{
    int m = degree + 1;
    vector<vector<double>> A(m, vector<double>(m));
    vector<double> B(m);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j)
            A[i][j] = simps_int([i, j](double x) { return phi(i, x) * phi(j, x); },
                a, b, n_int);
        B[i] = simps_int([&f, i](double x) { return f(x) * phi(i, x); },
            a, b, n_int);
    }
    return gauss_solve(A, B);
}

// ── Evaluate approximating polynomial ─────────────────────────────────────
double eval_approx(const vector<double>& coeffs, double x)
{
    double result = 0.0;
    for (int i = 0; i < (int)coeffs.size(); ++i)
        result += coeffs[i] * phi(i, x);
    return result;
}

// ── MSE (integrated squared error) ────────────────────────────────────────
double approx_mse(function<double(double)> f, const vector<double>& coeffs,
    double a, double b, int steps)
{
    double dx = (b - a) / steps;
    double Hn = 0.0;
    for (int s = 0; s <= steps; ++s) {
        double x = a + s * dx;
        double d = std::fabs(eval_approx(coeffs, x) - f(x));
        double w = (s == 0 || s == steps) ? 1.0 : (s % 2 == 0 ? 2.0 : 4.0);
        Hn += w * d * d;
    }
    return (dx / 3.0) * Hn;
}

// ── Max absolute error ─────────────────────────────────────────────────────
double approx_max_error(function<double(double)> f, const vector<double>& coeffs,
    double a, double b, int steps)
{
    double maxErr = 0.0;
    double dx = (b - a) / steps;
    for (int s = 0; s <= steps; ++s) {
        double x = a + s * dx;
        double d = std::fabs(eval_approx(coeffs, x) - f(x));
        if (d > maxErr) maxErr = d;
    }
    return maxErr;
}