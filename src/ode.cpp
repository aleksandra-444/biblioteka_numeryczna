#include "ode.h"
#include <cmath>
#include <stdexcept>
using std::function;

// ── Pojedyncze kroki ──────────────────────────────────────────────────────
double euler_step(RHS f, double y, double t, double h)
{
    return y + h * f(y, t);
}

double heun_step(RHS f, double y, double t, double h)
{
    double k1 = f(y, t);
    double tilde = y + h * k1;
    return y + (h / 2.0) * (k1 + f(tilde, t + h));
}

double midpoint_step(RHS f, double y, double t, double h)
{
    double mid = y + (h / 2.0) * f(y, t);
    return y + h * f(mid, t + h / 2.0);
}

double rk4_step(RHS f, double y, double t, double h)
{
    double k1 = h * f(y, t);
    double k2 = h * f(y + k1 / 2.0, t + h / 2.0);
    double k3 = h * f(y + k2 / 2.0, t + h / 2.0);
    double k4 = h * f(y + k3, t + h);
    return y + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// ── Pełne solvery ─────────────────────────────────────────────────────────
static double solve_generic(RHS f, double y0, double t_end, int N,
    int method)
{
    double h = t_end / N;
    double y = y0, t = 0.0;
    for (int i = 0; i < N; ++i) {
        switch (method) {
        case 0: y = euler_step(f, y, t, h); break;
        case 1: y = heun_step(f, y, t, h); break;
        case 2: y = midpoint_step(f, y, t, h); break;
        case 3: y = rk4_step(f, y, t, h); break;
        default: throw std::invalid_argument("unknown ODE method");
        }
        t += h;
    }
    return y;
}

double euler_solve(RHS f, double y0, double t_end, int N) { return solve_generic(f, y0, t_end, N, 0); }
double heun_solve(RHS f, double y0, double t_end, int N) { return solve_generic(f, y0, t_end, N, 1); }
double midpoint_solve(RHS f, double y0, double t_end, int N) { return solve_generic(f, y0, t_end, N, 2); }
double rk4_solve(RHS f, double y0, double t_end, int N) { return solve_generic(f, y0, t_end, N, 3); }

// ── MSE wzdłuż trajektorii ────────────────────────────────────────────────
double ode_mse(RHS f, function<double(double)> exact,
    double y0, double t_end, int N, int method)
{
    double h = t_end / N;
    double y = y0, t = 0.0;
    double sum = 0.0;
    for (int i = 0; i <= N; ++i) {
        double e = y - exact(t);
        sum += e * e;
        if (i < N) {
            switch (method) {
            case 0: y = euler_step(f, y, t, h); break;
            case 1: y = heun_step(f, y, t, h); break;
            case 2: y = midpoint_step(f, y, t, h); break;
            case 3: y = rk4_step(f, y, t, h); break;
            }
            t += h;
        }
    }
    return sum / (N + 1);
}