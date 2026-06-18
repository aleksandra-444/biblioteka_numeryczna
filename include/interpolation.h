#pragma once
/**
 * @file interpolation.h
 * @brief Interpolacja wielomianowa: metoda Lagrange'a i Newtona (różnice dzielone)
 *
 * Przykład:
 *   vector<double> xs = {0, 1, 2};
 *   vector<double> ys = {1, 3, 7};
 *   double val = lagrange_interpolate(1.5, xs, ys);
 */

#include <vector>
using std::vector;

/**
 * @brief Oblicza wartość wielomianu interpolacyjnego Lagrange'a w punkcie x
 * @param x  Punkt, w którym obliczana jest wartość
 * @param xs Węzły interpolacji (współrzędne x)
 * @param ys Wartości funkcji w węzłach (współrzędne y)
 * @return Wartość interpolowana L(x)
 */
double lagrange_interpolate(double x,
    const vector<double>& xs,
    const vector<double>& ys);

/**
 * @brief Buduje tablicę różnic dzielonych Newtona
 * @param xs Węzły interpolacji
 * @param ys Wartości funkcji w węzłach
 * @return Dwuwymiarowa tablica f[i][j] różnic dzielonych
 */
vector<vector<double>> newton_divided_differences(const vector<double>& xs,
    const vector<double>& ys);

/**
 * @brief Oblicza wartość wielomianu interpolacyjnego Newtona przy użyciu gotowej tablicy
 * @param x  Punkt, w którym obliczana jest wartość
 * @param xs Węzły interpolacji
 * @param dd Tablica różnic dzielonych z newton_divided_differences()
 * @return Wartość interpolowana N(x)
 */
double newton_interpolate(double x,
    const vector<double>& xs,
    const vector<vector<double>>& dd);

/**
 * @brief Oblicza średni błąd kwadratowy interpolacji dla wszystkich punktów danych
 * @param xs_all Wszystkie punkty x
 * @param ys_all Odpowiadające wartości funkcji
 * @param xs_nodes Wybrane węzły interpolacji (podzbiór)
 * @param ys_nodes Wartości funkcji w węzłach
 * @return Wartość MSE
 */
double interpolation_mse(const vector<double>& xs_all,
    const vector<double>& ys_all,
    const vector<double>& xs_nodes,
    const vector<double>& ys_nodes);

/**
 * @brief Oblicza wartość wielomianu w postaci naturalnej schematem Hornera
 * @param a  Współczynniki a[0] + a[1]*x + ... + a[n]*x^n
 * @param x  Punkt obliczenia
 * @return   Wartość wielomianu
 *
 * Przykład:
 *   vector<double> a = {1, 0, 1}; // 1 + x^2
 *   double v = horner_eval(a, 3.0); // zwraca 10.0
 */
double horner_eval(const vector<double>& a, double x);