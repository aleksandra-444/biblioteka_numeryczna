#pragma once
/**
 * @file approximation.h
 * @brief Least-squares polynomial approximation using power basis phi_i(x) = x^i
 *
 * Builds and solves the normal equations via Gaussian elimination.
 *
 * Example:
 *   // Approximate f(x) = (x-1)*exp(sin(x^2)) on [-2.5, 1.0] with degree 10
 *   auto f = [](double x){ return (x-1)*exp(sin(x*x)); };
 *   auto coeffs = least_squares_poly(f, -2.5, 1.0, 10);
 *   double val = eval_approx(coeffs, 0.5);
 */

#include <vector>
#include <functional>
using std::vector;
using std::function;

/**
 * @brief Compute least-squares polynomial approximation coefficients
 * @param f      Target function to approximate
 * @param a      Left endpoint of interval
 * @param b      Right endpoint of interval
 * @param degree Degree of approximating polynomial (basis size = degree+1)
 * @param n_int  Number of integration sub-intervals for inner products (default 1000)
 * @return Coefficient vector c such that P(x) = c[0] + c[1]*x + ... + c[degree]*x^degree
 */
vector<double> least_squares_poly(function<double(double)> f,
    double a, double b,
    int degree,
    int n_int = 1000);

/**
 * @brief Evaluate approximating polynomial at x
 * @param coeffs Coefficient vector from least_squares_poly()
 * @param x      Evaluation point
 * @return       P(x)
 */
double eval_approx(const vector<double>& coeffs, double x);

/**
 * @brief Compute mean-squared approximation error (Hn) over [a,b]
 * @param f      Original function
 * @param coeffs Approximation coefficients
 * @param a      Left endpoint
 * @param b      Right endpoint
 * @param steps  Number of evaluation steps (default 1000)
 * @return Integrated squared error
 */
double approx_mse(function<double(double)> f,
    const vector<double>& coeffs,
    double a, double b,
    int steps = 1000);

/**
 * @brief Compute maximum absolute error of approximation over [a,b]
 * @param f      Original function
 * @param coeffs Approximation coefficients
 * @param a      Left endpoint
 * @param b      Right endpoint
 * @param steps  Number of evaluation steps
 * @return Maximum |f(x) - P(x)|
 */
double approx_max_error(function<double(double)> f,
    const vector<double>& coeffs,
    double a, double b,
    int steps = 1000);