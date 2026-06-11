#pragma once
/**
 * @file nonlinear.h
 * @brief Root-finding methods: Bisection, Regula Falsi, Newton, Secant
 *
 * Example:
 *   auto f  = [](double x){ return x*x - 2.0; };
 *   auto df = [](double x){ return 2.0*x; };
 *   double root = newton_method(f, df, 2.0); // approx sqrt(2)
 */

#include <functional>
using std::function;

/** Maximum number of iterations used by all methods */
static const int NL_MAX_ITER = 150;
/** Default convergence tolerance */
static const double NL_EPS = 1e-9;

/**
 * @brief Bisection method
 * @param f  Continuous function with sign change on [a, b]
 * @param a  Left bracket
 * @param b  Right bracket
 * @param eps Convergence tolerance (default NL_EPS)
 * @return Approximate root
 */
double bisection(function<double(double)> f, double a, double b,
    double eps = NL_EPS);

/**
 * @brief Regula Falsi (False Position) method
 * @param f  Continuous function with sign change on [a, b]
 * @param a  Left bracket
 * @param b  Right bracket
 * @param eps Convergence tolerance
 * @return Approximate root
 */
double regula_falsi(function<double(double)> f, double a, double b,
    double eps = NL_EPS);

/**
 * @brief Newton-Raphson method with analytical derivative
 * @param f   Function
 * @param df  Derivative of f
 * @param x0  Initial guess
 * @param eps Convergence tolerance
 * @return Approximate root
 */
double newton_method(function<double(double)> f,
    function<double(double)> df,
    double x0,
    double eps = NL_EPS);

/**
 * @brief Newton-Raphson method with numerical derivative (central differences)
 * @param f   Function
 * @param x0  Initial guess
 * @param eps Convergence tolerance
 * @param h   Step for numerical differentiation (default 1e-5)
 * @return Approximate root
 */
double newton_numeric(function<double(double)> f, double x0,
    double eps = NL_EPS, double h = 1e-5);

/**
 * @brief Secant method
 * @param f   Function
 * @param x0  First initial guess
 * @param x1  Second initial guess
 * @param eps Convergence tolerance
 * @return Approximate root
 */
double secant_method(function<double(double)> f, double x0, double x1,
    double eps = NL_EPS);

/**
 * @brief Count iterations used by bisection to reach tolerance eps
 */
int bisection_iterations(function<double(double)> f, double a, double b,
    double eps = NL_EPS);

/**
 * @brief Count iterations used by Newton method to reach tolerance eps
 */
int newton_iterations(function<double(double)> f,
    function<double(double)> df,
    double x0, double eps = NL_EPS);