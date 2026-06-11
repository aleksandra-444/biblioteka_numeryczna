#pragma once
/**
 * @file linear_systems.h
 * @brief Solving systems of linear equations: LU decomposition, Gaussian elimination
 *
 * Provides LU decomposition (Doolittle), forward/backward substitution,
 * and Gaussian elimination with partial pivoting.
 *
 * Example:
 *   vector<vector<double>> A = {{2,1},{5,3}};
 *   vector<double> b = {8, 19};
 *   auto x = gauss_solve(A, b);
 *   // x = {5, -2}
 */

#include <vector>
#include <string>
using std::vector;
using std::string;

/**
 * @brief LU decomposition (Doolittle method, no pivoting)
 * @param A Input square matrix (n x n)
 * @param L Output lower triangular matrix (L[i][i] = 1)
 * @param U Output upper triangular matrix
 * @param n Matrix size
 */
void lu_decompose(const vector<vector<double>>& A,
    vector<vector<double>>& L,
    vector<vector<double>>& U,
    int n);

/**
 * @brief Forward substitution: solve Lz = b
 * @param L Lower triangular matrix
 * @param b Right-hand side vector
 * @param n System size
 * @return Solution vector z
 */
vector<double> forward_substitution(const vector<vector<double>>& L,
    const vector<double>& b,
    int n);

/**
 * @brief Backward substitution: solve Ux = z
 * @param U Upper triangular matrix
 * @param z Right-hand side vector
 * @param n System size
 * @return Solution vector x
 */
vector<double> backward_substitution(const vector<vector<double>>& U,
    const vector<double>& z,
    int n);

/**
 * @brief Solve Ax = b using LU decomposition
 * @param A Square matrix (n x n)
 * @param b Right-hand side vector
 * @return Solution vector x
 */
vector<double> lu_solve(const vector<vector<double>>& A,
    const vector<double>& b);

/**
 * @brief Gaussian elimination with partial pivoting
 * @param A Square matrix (n x n), will NOT be modified (copy used internally)
 * @param b Right-hand side vector
 * @return Solution vector x
 */
vector<double> gauss_solve(vector<vector<double>> A, vector<double> b);

/**
 * @brief Analyse a linear system for consistency and uniqueness
 * @param A Coefficient matrix
 * @param b Right-hand side vector
 * @param n System size
 * @return "OZNACZONY", "NIEOZNACZONY", or "SPRZECZNY"
 */
string analyze_system(vector<vector<double>> A, vector<double> b, int n);

/**
 * @brief Compute residual norm ||Ax - b||_2
 * @param A Coefficient matrix
 * @param x Solution vector
 * @param b Right-hand side vector
 * @param n System size
 * @return Euclidean norm of the residual
 */
double residual_norm(const vector<vector<double>>& A,
    const vector<double>& x,
    const vector<double>& b,
    int n);