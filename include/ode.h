#pragma once
/**
 * @file ode.h
 * @brief Numeryczne metody całkowania równań różniczkowych: Euler, Heun, Punkt środkowy, RK4
 *
 * Wszystkie solvery całkują skalarne RRZ  dy/dt = f(y, t)  od t=0 do t=t_end
 * używając N równych kroków i zwracają wartość końcową.
 *
 * Przykład:
 *   // Prawo stygnięcia: dT/dt = -alpha*T^4
 *   auto f = [](double T, double t){ return -3e-12 * (T * T * T * T); };
 *   double T_end = euler_solve(f, 5903.0, 5903.0, 10000);
 */

#include <functional>
#include <vector>
     using std::function;
 using std::vector;

 /** Typ prawej strony równania: f(y, t) */
 using RHS = function<double(double, double)>;

 /**
  * @brief Pojedynczy krok metody Eulera
  * @param f  Funkcja prawej strony f(y, t)
  * @param y  Aktualna wartość
  * @param t  Aktualny czas
  * @param h  Rozmiar kroku
  * @return   y w chwili t+h
  */
 double euler_step(RHS f, double y, double t, double h);

 /**
  * @brief Pojedynczy krok metody Heuna (ulepszony Euler / RK2)
  */
 double heun_step(RHS f, double y, double t, double h);

 /**
  * @brief Pojedynczy krok metody punktu środkowego
  */
 double midpoint_step(RHS f, double y, double t, double h);

 /**
  * @brief Pojedynczy krok klasycznej metody RK4
  */
 double rk4_step(RHS f, double y, double t, double h);

 /**
  * @brief Całkuje RRZ metodą Eulera przez N kroków
  * @param f      Funkcja prawej strony
  * @param y0     Wartość początkowa
  * @param t_end  Czas końcowy (czas startowy = 0)
  * @param N      Liczba kroków
  * @return Wartość końcowa y(t_end)
  */
 double euler_solve(RHS f, double y0, double t_end, int N);

 /** @brief Całkuje RRZ metodą Heuna */
 double heun_solve(RHS f, double y0, double t_end, int N);

 /** @brief Całkuje RRZ metodą punktu środkowego */
 double midpoint_solve(RHS f, double y0, double t_end, int N);

 /** @brief Całkuje RRZ metodą RK4 */
 double rk4_solve(RHS f, double y0, double t_end, int N);

 /**
  * @brief Oblicza MSE między rozwiązaniem numerycznym a dokładnym w N+1 punktach
  * @param f       Prawa strona RRZ
  * @param exact   Rozwiązanie dokładne exact(t)
  * @param y0      Wartość początkowa
  * @param t_end   Czas końcowy
  * @param N       Liczba kroków
  * @param method  0=Euler, 1=Heun, 2=Punkt środkowy, 3=RK4
  * @return MSE
  */
 double ode_mse(RHS f,
     function<double(double)> exact,
     double y0, double t_end, int N,
     int method);