#include "integration.h"
#include <vector>
#include <stdexcept>
using std::vector;
using std::function;

// ── Trapezoidal rule ───────────────────────────────────────────────────────
double trapezoidal(function<double(double)> f, double a, double b, int m)
{
    double h = (b - a) / m;
    double s = 0.5 * (f(a) + f(b));
    for (int i = 1; i < m; ++i) s += f(a + i * h);
    return s * h;
}

// ── Simpson's rule ─────────────────────────────────────────────────────────
double simpson(function<double(double)> f, double a, double b, int m)
{
    if (m % 2 != 0) ++m; // enforce even
    double h = (b - a) / m;
    double s = f(a) + f(b);
    for (int i = 1; i < m; ++i)
        s += (i % 2 == 0 ? 2.0 : 4.0) * f(a + i * h);
    return (h / 3.0) * s;
}

// ── Gauss-Legendre nodes & weights (reference interval [-1,1]) ─────────────
static void gauss_nodes_weights(int n, vector<double>& xi, vector<double>& wi)
{
    if (n == 2) {
        xi = { -0.5773502691896257,  0.5773502691896257 };
        wi = { 1.0,                  1.0 };
    }
    else if (n == 3) {
        xi = { 0.0,                 -0.7745966692414834, 0.7745966692414834 };
        wi = { 8.0 / 9.0,              5.0 / 9.0,            5.0 / 9.0 };
    }
    else if (n == 4) {
        xi = { -0.3399810435848563,   0.3399810435848563,
              -0.8611363115940526,   0.8611363115940526 };
        wi = { 0.6521451548625461,   0.6521451548625461,
               0.3478548451374538,   0.3478548451374538 };
    }
    else if (n == 5) {
        xi = { 0.0,  0.538469, -0.538469,  0.90618, -0.90618 };
        wi = { 0.56889, 0.478629, 0.478629, 0.236927, 0.236927 };
    }
    else {
        throw std::invalid_argument("gauss_legendre: n_nodes must be 2, 3, 4, or 5");
    }
}

// ── Single-interval Gauss-Legendre ────────────────────────────────────────
double gauss_legendre(function<double(double)> f, double a, double b, int n_nodes)
{
    vector<double> xi, wi;
    gauss_nodes_weights(n_nodes, xi, wi);
    double mid = 0.5 * (a + b);
    double half = 0.5 * (b - a);
    double s = 0.0;
    for (int i = 0; i < n_nodes; ++i)
        s += wi[i] * f(half * xi[i] + mid);
    return half * s;
}

// ── Composite Gauss-Legendre ──────────────────────────────────────────────
double gauss_legendre_composite(function<double(double)> f,
    double a, double b, int m, int n_nodes)
{
    double h = (b - a) / m;
    double s = 0.0;
    for (int i = 0; i < m; ++i)
        s += gauss_legendre(f, a + i * h, a + (i + 1) * h, n_nodes);
    return s;
}