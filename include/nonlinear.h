#pragma once
/**
 * @file nonlinear.h
 * @brief Metody znajdowania pierwiastków: Bisekcja, Regula Falsi, Newton, Sieczna
 *
 * Przykład:
 *   auto f  = [](double x){ return x*x - 2.0; };
 *   auto df = [](double x){ return 2.0*x; };
 *   double root = newton_method(f, df, 2.0); // przybliżenie sqrt(2)
 */

#include <functional>
using std::function;

/** Maksymalna liczba iteracji używana przez wszystkie metody */
static const int NL_MAX_ITER = 150;
/** Domyślna tolerancja zbieżności */
static const double NL_EPS = 1e-9;

/**
 * @brief Metoda bisekcji
 * @param f  Ciągła funkcja ze zmianą znaku na [a, b]
 * @param a  Lewa granica przedziału
 * @param b  Prawa granica przedziału
 * @param eps Tolerancja zbieżności (domyślnie NL_EPS)
 * @return Przybliżony pierwiastek
 */
double bisection(function<double(double)> f, double a, double b,
    double eps = NL_EPS);

/**
 * @brief Metoda Regula Falsi (fałszywego położenia)
 * @param f  Ciągła funkcja ze zmianą znaku na [a, b]
 * @param a  Lewa granica przedziału
 * @param b  Prawa granica przedziału
 * @param eps Tolerancja zbieżności
 * @return Przybliżony pierwiastek
 */
double regula_falsi(function<double(double)> f, double a, double b,
    double eps = NL_EPS);

/**
 * @brief Metoda Newtona-Raphsona z analityczną pochodną
 * @param f   Funkcja
 * @param df  Pochodna funkcji f
 * @param x0  Punkt startowy
 * @param eps Tolerancja zbieżności
 * @return Przybliżony pierwiastek
 */
double newton_method(function<double(double)> f,
    function<double(double)> df,
    double x0,
    double eps = NL_EPS);

/**
 * @brief Metoda Newtona-Raphsona z numeryczną pochodną (różnice centralne)
 * @param f   Funkcja
 * @param x0  Punkt startowy
 * @param eps Tolerancja zbieżności
 * @param h   Krok różniczkowania numerycznego (domyślnie 1e-5)
 * @return Przybliżony pierwiastek
 */
double newton_numeric(function<double(double)> f, double x0,
    double eps = NL_EPS, double h = 1e-5);

/**
 * @brief Metoda siecznych
 * @param f   Funkcja
 * @param x0  Pierwsze przybliżenie
 * @param x1  Drugie przybliżenie
 * @param eps Tolerancja zbieżności
 * @return Przybliżony pierwiastek
 */
double secant_method(function<double(double)> f, double x0, double x1,
    double eps = NL_EPS);

/**
 * @brief Zlicza iteracje bisekcji potrzebne do osiągnięcia tolerancji eps
 */
int bisection_iterations(function<double(double)> f, double a, double b,
    double eps = NL_EPS);

/**
 * @brief Zlicza iteracje metody Newtona potrzebne do osiągnięcia tolerancji eps
 */
int newton_iterations(function<double(double)> f,
    function<double(double)> df,
    double x0, double eps = NL_EPS);