#pragma once
/**
 * @file integration.h
 * @brief Numerical integration: Trapezoidal, Simpson, Gauss-Legendre quadrature
 *
 * Example:
 *   auto f = [](double x){ return x*x; };
 *   double I = simpson(f, 0.0, 1.0, 1000); // approx 1/3
 */

#include <functional>
using std::function;

/**
 * @brief Composite trapezoidal rule
 * @param f  Integrand function
 * @param a  Left endpoint
 * @param b  Right endpoint
 * @param m  Number of sub-intervals (must be >= 1)
 * @return   Approximate integral
 */
double trapezoidal(function<double(double)> f, double a, double b, int m);

/**
 * @brief Composite Simpson's rule (requires even m)
 * @param f  Integrand function
 * @param a  Left endpoint
 * @param b  Right endpoint
 * @param m  Number of sub-intervals (must be even and >= 2)
 * @return   Approximate integral
 */
double simpson(function<double(double)> f, double a, double b, int m);

/**
 * @brief Single-interval Gauss-Legendre quadrature on [a, b]
 * @param f       Integrand function
 * @param a       Left endpoint
 * @param b       Right endpoint
 * @param n_nodes Number of Gauss nodes (2, 3, 4, or 5)
 * @return        Approximate integral
 */
double gauss_legendre(function<double(double)> f, double a, double b, int n_nodes);

/**
 * @brief Composite Gauss-Legendre quadrature: subdivide [a,b] into m intervals
 * @param f       Integrand function
 * @param a       Left endpoint
 * @param b       Right endpoint
 * @param m       Number of sub-intervals
 * @param n_nodes Number of Gauss nodes per interval (2, 3, 4, or 5)
 * @return        Approximate integral
 */
double gauss_legendre_composite(function<double(double)> f,
    double a, double b,
    int m, int n_nodes);