#pragma once
/**
 * @file approximation.h
 * @brief Aproksymacja wielomianowa metodą najmniejszych kwadratów, baza potęgowa phi_i(x) = x^i
 *
 * Buduje i rozwiązuje równania normalne za pomocą eliminacji Gaussa.
 *
 * Przykład:
 *   // Aproksymacja f(x) = (x-1)*exp(sin(x^2)) na [-2.5, 1.0] wielomianem stopnia 10
 *   auto f = [](double x){ return (x-1)*exp(sin(x*x)); };
 *   auto coeffs = least_squares_poly(f, -2.5, 1.0, 10);
 *   double val = eval_approx(coeffs, 0.5);
 */

#include <vector>
#include <functional>
using std::vector;
using std::function;

/**
 * @brief Oblicza współczynniki aproksymacji wielomianowej metodą najmniejszych kwadratów
 * @param f      Aproksymowana funkcja
 * @param a      Lewy koniec przedziału
 * @param b      Prawy koniec przedziału
 * @param degree Stopień wielomianu aproksymującego (rozmiar bazy = degree+1)
 * @param n_int  Liczba podprzedziałów całkowania dla iloczynów skalarnych (domyślnie 1000)
 * @return Wektor współczynników c, gdzie P(x) = c[0] + c[1]*x + ... + c[degree]*x^degree
 */
vector<double> least_squares_poly(function<double(double)> f,
    double a, double b,
    int degree,
    int n_int = 1000);

/**
 * @brief Oblicza wartość wielomianu aproksymującego w punkcie x
 * @param coeffs Wektor współczynników z least_squares_poly()
 * @param x      Punkt obliczenia
 * @return       P(x)
 */
double eval_approx(const vector<double>& coeffs, double x);

/**
 * @brief Oblicza średni błąd kwadratowy aproksymacji (Hn) na przedziale [a,b]
 * @param f      Oryginalna funkcja
 * @param coeffs Współczynniki aproksymacji
 * @param a      Lewy koniec przedziału
 * @param b      Prawy koniec przedziału
 * @param steps  Liczba kroków obliczeniowych (domyślnie 1000)
 * @return Scałkowany błąd kwadratowy
 */
double approx_mse(function<double(double)> f,
    const vector<double>& coeffs,
    double a, double b,
    int steps = 1000);

/**
 * @brief Oblicza maksymalny błąd bezwzględny aproksymacji na przedziale [a,b]
 * @param f      Oryginalna funkcja
 * @param coeffs Współczynniki aproksymacji
 * @param a      Lewy koniec przedziału
 * @param b      Prawy koniec przedziału
 * @param steps  Liczba kroków obliczeniowych
 * @return Maksimum |f(x) - P(x)|
 */
double approx_max_error(function<double(double)> f,
    const vector<double>& coeffs,
    double a, double b,
    int steps = 1000);