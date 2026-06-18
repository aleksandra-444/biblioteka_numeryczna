#include "interpolation.h"
#include <cmath>
using std::vector;

// ── Interpolacja Lagrange'a ────────────────────────────────────────────────
double lagrange_interpolate(double x, const vector<double>& xs, const vector<double>& ys)
{
    int n = (int)xs.size();
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        double L = 1.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) L *= (x - xs[j]) / (xs[i] - xs[j]);
        }
        result += ys[i] * L;
    }
    return result;
}

// ── Różnice dzielone Newtona ───────────────────────────────────────────────
vector<vector<double>> newton_divided_differences(const vector<double>& xs,
    const vector<double>& ys)
{
    int n = (int)xs.size();
    vector<vector<double>> f(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) f[i][0] = ys[i];
    for (int j = 1; j < n; ++j)
        for (int i = 0; i < n - j; ++i)
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (xs[i + j] - xs[i]);
    return f;
}

// ── Obliczanie wartości wielomianu interpolacyjnego Newtona ───────────────
double newton_interpolate(double x, const vector<double>& xs,
    const vector<vector<double>>& dd)
{
    int n = (int)xs.size();
    double result = dd[0][n - 1];
    for (int i = n - 2; i >= 0; --i)
        result = dd[0][i] + (x - xs[i]) * result;
    return result;
}

// ── Średni błąd kwadratowy interpolacji ───────────────────────────────────
double interpolation_mse(const vector<double>& xs_all, const vector<double>& ys_all,
    const vector<double>& xs_nodes, const vector<double>& ys_nodes)
{
    auto dd = newton_divided_differences(xs_nodes, ys_nodes);
    double sum = 0.0;
    int N = (int)xs_all.size();
    for (int i = 0; i < N; ++i) {
        double val = newton_interpolate(xs_all[i], xs_nodes, dd);
        double err = val - ys_all[i];
        sum += err * err;
    }
    return sum / N;
}

// ── Schemat Hornera ────────────────────────────────────────────────────────
double horner_eval(const vector<double>& a, double x)
{
    int n = (int)a.size() - 1;
    double result = a[n];
    for (int i = n - 1; i >= 0; --i)
        result = result * x + a[i];
    return result;
}