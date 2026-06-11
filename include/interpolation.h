#pragma once
/**
 * @file interpolation.h
 * @brief Polynomial interpolation: Lagrange and Newton (divided differences)
 *
 * Example:
 *   vector<double> xs = {0, 1, 2};
 *   vector<double> ys = {1, 3, 7};
 *   double val = lagrange_interpolate(1.5, xs, ys);
 */

#include <vector>
using std::vector;

/**
 * @brief Evaluate Lagrange interpolating polynomial at point x
 * @param x  Point at which to evaluate
 * @param xs Interpolation nodes (x-coordinates)
 * @param ys Function values at nodes (y-coordinates)
 * @return Interpolated value L(x)
 */
double lagrange_interpolate(double x,
    const vector<double>& xs,
    const vector<double>& ys);

/**
 * @brief Build Newton divided-difference table
 * @param xs Interpolation nodes
 * @param ys Function values at nodes
 * @return 2D table f[i][j] of divided differences
 */
vector<vector<double>> newton_divided_differences(const vector<double>& xs,
    const vector<double>& ys);

/**
 * @brief Evaluate Newton interpolating polynomial using precomputed table
 * @param x  Point at which to evaluate
 * @param xs Interpolation nodes
 * @param dd Divided-difference table from newton_divided_differences()
 * @return Interpolated value N(x)
 */
double newton_interpolate(double x,
    const vector<double>& xs,
    const vector<vector<double>>& dd);

/**
 * @brief Compute mean-squared error of interpolation over all given data points
 * @param xs_all All x data points
 * @param ys_all Corresponding function values
 * @param xs_nodes Chosen interpolation nodes (subset)
 * @param ys_nodes Function values at nodes
 * @return MSE value
 */
double interpolation_mse(const vector<double>& xs_all,
    const vector<double>& ys_all,
    const vector<double>& xs_nodes,
    const vector<double>& ys_nodes);

/**
 * @brief Evaluate polynomial in natural form using Horner's scheme
 * @param a  Coefficients a[0] + a[1]*x + ... + a[n]*x^n
 * @param x  Evaluation point
 * @return   Polynomial value
 *
 * Example:
 *   vector<double> a = {1, 0, 1}; // 1 + x^2
 *   double v = horner_eval(a, 3.0); // returns 10.0
 */
double horner_eval(const vector<double>& a, double x);