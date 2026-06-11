#include "linear_systems.h"
#include <cmath>
#include <stdexcept>
using std::vector;
using std::string;
using std::swap;
using std::fabs;

// ── LU decomposition (Doolittle, no pivoting) ─────────────────────────────
void lu_decompose(const vector<vector<double>>& A,
    vector<vector<double>>& L,
    vector<vector<double>>& U,
    int n)
{
    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) L[i][i] = 1.0;

    for (int k = 0; k < n; ++k) {
        // Fill row k of U
        for (int j = k; j < n; ++j) {
            double s = 0.0;
            for (int p = 0; p < k; ++p) s += L[k][p] * U[p][j];
            U[k][j] = A[k][j] - s;
        }
        if (fabs(U[k][k]) < 1e-12) U[k][k] = 1e-12; // guard

        // Fill column k of L
        for (int i = k + 1; i < n; ++i) {
            double s = 0.0;
            for (int p = 0; p < k; ++p) s += L[i][p] * U[p][k];
            L[i][k] = (A[i][k] - s) / U[k][k];
        }
    }
}

// ── Forward substitution  Lz = b ──────────────────────────────────────────
vector<double> forward_substitution(const vector<vector<double>>& L,
    const vector<double>& b, int n)
{
    vector<double> z(n);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < i; ++j) s += L[i][j] * z[j];
        z[i] = b[i] - s;
    }
    return z;
}

// ── Backward substitution  Ux = z ─────────────────────────────────────────
vector<double> backward_substitution(const vector<vector<double>>& U,
    const vector<double>& z, int n)
{
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double s = 0.0;
        for (int j = i + 1; j < n; ++j) s += U[i][j] * x[j];
        if (fabs(U[i][i]) < 1e-15)
            throw std::runtime_error("lu_solve: singular matrix");
        x[i] = (z[i] - s) / U[i][i];
    }
    return x;
}

// ── Solve via LU ──────────────────────────────────────────────────────────
vector<double> lu_solve(const vector<vector<double>>& A, const vector<double>& b)
{
    int n = (int)b.size();
    vector<vector<double>> L, U;
    lu_decompose(A, L, U, n);
    auto z = forward_substitution(L, b, n);
    return backward_substitution(U, z, n);
}

// ── Gaussian elimination with partial pivoting ────────────────────────────
vector<double> gauss_solve(vector<vector<double>> A, vector<double> b)
{
    int n = (int)b.size();
    for (int i = 0; i < n; ++i) {
        // Find pivot
        int maxRow = i;
        for (int k = i + 1; k < n; ++k)
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) maxRow = k;
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);
        if (fabs(A[i][i]) < 1e-12) continue;
        // Eliminate
        for (int k = i + 1; k < n; ++k) {
            double c = -A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) A[k][j] += c * A[i][j];
            b[k] += c * b[i];
        }
    }
    // Back-substitute
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        if (fabs(A[i][i]) < 1e-12) continue;
        x[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; --k) b[k] -= A[k][i] * x[i];
    }
    return x;
}

// ── System analysis ───────────────────────────────────────────────────────
string analyze_system(vector<vector<double>> A, vector<double> b, int n)
{
    // Build augmented matrix
    vector<vector<double>> M(n, vector<double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }
    int row = 0;
    for (int col = 0; col < n && row < n; ++col) {
        int pivot = row;
        for (int i = row; i < n; ++i)
            if (fabs(M[i][col]) > fabs(M[pivot][col])) pivot = i;
        if (fabs(M[pivot][col]) < 1e-12) continue;
        swap(M[pivot], M[row]);
        for (int i = row + 1; i < n; ++i) {
            double f = M[i][col] / M[row][col];
            for (int j = col; j <= n; ++j) M[i][j] -= f * M[row][j];
        }
        ++row;
    }
    int rankA = 0, rankAB = 0;
    for (int i = 0; i < n; ++i) {
        bool zA = true, zAB = true;
        for (int j = 0; j < n; ++j)  if (fabs(M[i][j]) > 1e-10) zA = false;
        for (int j = 0; j <= n; ++j) if (fabs(M[i][j]) > 1e-10) zAB = false;
        if (!zA)  ++rankA;
        if (!zAB) ++rankAB;
    }
    if (rankA != rankAB) return "SPRZECZNY";
    if (rankA < n)       return "NIEOZNACZONY";
    return "OZNACZONY";
}

// ── Residual norm ─────────────────────────────────────────────────────────
double residual_norm(const vector<vector<double>>& A,
    const vector<double>& x,
    const vector<double>& b, int n)
{
    double err = 0.0;
    for (int i = 0; i < n; ++i) {
        double row = 0.0;
        for (int j = 0; j < n; ++j) row += A[i][j] * x[j];
        err += (row - b[i]) * (row - b[i]);
    }
    return std::sqrt(err);
}