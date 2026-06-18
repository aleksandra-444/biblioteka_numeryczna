#pragma once
/**
 * @file integration.h
 * @brief Całkowanie numeryczne: reguła trapezów, Simpsona, kwadratura Gaussa-Legendre'a
 *
 * Przykład:
 *   auto f = [](double x){ return x*x; };
 *   double I = simpson(f, 0.0, 1.0, 1000); // przybliżenie 1/3
 */

#include <functional>
using std::function;

/**
 * @brief Złożona reguła trapezów
 * @param f  Funkcja podcałkowa
 * @param a  Lewy koniec przedziału
 * @param b  Prawy koniec przedziału
 * @param m  Liczba podprzedziałów (>= 1)
 * @return   Przybliżona wartość całki
 */
double trapezoidal(function<double(double)> f, double a, double b, int m);

/**
 * @brief Złożona reguła Simpsona (wymaga parzystego m)
 * @param f  Funkcja podcałkowa
 * @param a  Lewy koniec przedziału
 * @param b  Prawy koniec przedziału
 * @param m  Liczba podprzedziałów (parzysta, >= 2)
 * @return   Przybliżona wartość całki
 */
double simpson(function<double(double)> f, double a, double b, int m);

/**
 * @brief Kwadratura Gaussa-Legendre'a na pojedynczym przedziale [a, b]
 * @param f       Funkcja podcałkowa
 * @param a       Lewy koniec przedziału
 * @param b       Prawy koniec przedziału
 * @param n_nodes Liczba węzłów Gaussa (2, 3, 4 lub 5)
 * @return        Przybliżona wartość całki
 */
double gauss_legendre(function<double(double)> f, double a, double b, int n_nodes);

/**
 * @brief Złożona kwadratura Gaussa-Legendre'a: podział [a,b] na m podprzedziałów
 * @param f       Funkcja podcałkowa
 * @param a       Lewy koniec przedziału
 * @param b       Prawy koniec przedziału
 * @param m       Liczba podprzedziałów
 * @param n_nodes Liczba węzłów Gaussa na podprzedzial (2, 3, 4 lub 5)
 * @return        Przybliżona wartość całki
 */
double gauss_legendre_composite(function<double(double)> f,
    double a, double b,
    int m, int n_nodes);