#include "nonlinear.h"
#include <cmath>
#include <stdexcept>

using std::function;

double bisection(function<double(double)> f, double a, double b, double eps)
{
    double fa = f(a), fb = f(b);
    if (fa * fb >= 0.0)
        throw std::invalid_argument("bisection: f(a) and f(b) must have opposite signs");

    for (int i = 0; i < NL_MAX_ITER; ++i) {
        double mid = 0.5 * (a + b);
        double fmid = f(mid);
        if (std::fabs(b - a) < eps || fmid == 0.0)
            return mid;
        if (fa * fmid < 0.0) { b = mid; fb = fmid; }
        else                  { a = mid; fa = fmid; }
    }
    return 0.5 * (a + b);
}

double regula_falsi(function<double(double)> f, double a, double b, double eps)
{
    double fa = f(a), fb = f(b);
    if (fa * fb >= 0.0)
        throw std::invalid_argument("regula_falsi: f(a) and f(b) must have opposite signs");

    double c = a;
    for (int i = 0; i < NL_MAX_ITER; ++i) {
        c = a - fa * (b - a) / (fb - fa);
        double fc = f(c);
        if (std::fabs(fc) < eps)
            return c;
        if (fa * fc < 0.0) { b = c; fb = fc; }
        else               { a = c; fa = fc; }
    }
    return c;
}

double newton_method(function<double(double)> f, function<double(double)> df,
                     double x0, double eps)
{
    double x = x0;
    for (int i = 0; i < NL_MAX_ITER; ++i) {
        double fx = f(x);
        double dfx = df(x);
        if (std::fabs(dfx) < 1e-15)
            break;
        double xn = x - fx / dfx;
        if (std::fabs(xn - x) < eps)
            return xn;
        x = xn;
    }
    return x;
}

double newton_numeric(function<double(double)> f, double x0, double eps, double h)
{
    double x = x0;
    for (int i = 0; i < NL_MAX_ITER; ++i) {
        double fx = f(x);
        double dfx = (f(x + h) - f(x - h)) / (2.0 * h);
        if (std::fabs(dfx) < 1e-15)
            break;
        double xn = x - fx / dfx;
        if (std::fabs(xn - x) < eps)
            return xn;
        x = xn;
    }
    return x;
}

double secant_method(function<double(double)> f, double x0, double x1, double eps)
{
    double fx0 = f(x0), fx1 = f(x1);
    for (int i = 0; i < NL_MAX_ITER; ++i) {
        double denom = fx1 - fx0;
        if (std::fabs(denom) < 1e-15)
            break;
        double x2 = x1 - fx1 * (x1 - x0) / denom;
        if (std::fabs(x2 - x1) < eps)
            return x2;
        x0 = x1; fx0 = fx1;
        x1 = x2; fx1 = f(x2);
    }
    return x1;
}

int bisection_iterations(function<double(double)> f, double a, double b, double eps)
{
    double fa = f(a), fb = f(b);
    if (fa * fb >= 0.0)
        throw std::invalid_argument("bisection_iterations: f(a) and f(b) must have opposite signs");

    int count = 0;
    for (int i = 0; i < NL_MAX_ITER; ++i) {
        double mid = 0.5 * (a + b);
        double fmid = f(mid);
        ++count;
        if (std::fabs(b - a) < eps || fmid == 0.0)
            return count;
        if (fa * fmid < 0.0) { b = mid; fb = fmid; }
        else                  { a = mid; fa = fmid; }
    }
    return count;
}

int newton_iterations(function<double(double)> f, function<double(double)> df,
                      double x0, double eps)
{
    double x = x0;
    for (int i = 0; i < NL_MAX_ITER; ++i) {
        double fx = f(x);
        double dfx = df(x);
        if (std::fabs(dfx) < 1e-15)
            return i;
        double xn = x - fx / dfx;
        if (std::fabs(xn - x) < eps)
            return i + 1;
        x = xn;
    }
    return NL_MAX_ITER;
}