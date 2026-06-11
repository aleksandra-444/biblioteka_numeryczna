#pragma once
/**
 * @file ode.h
 * @brief Numerical ODE solvers: Euler, Heun, Midpoint, Runge-Kutta 4
 *
 * All solvers advance a scalar ODE  dT/dt = f(T, t)  from t=0 to t=t_end
 * using N equal steps and return the final value.
 *
 * Example:
 *   // Cooling law: dT/dt = -alpha*(T^4 - beta)
 *   auto f = [](double T, double t){ return -3e-12 * (T * T * T * T); };
 *double T_end = euler_step_n(f, 5903.0, 5903.0, 10000);
 */

#include <functional>
#include <vector>
     using std::function;
 using std::vector;

 /** Right-hand side type: f(y, t) */
 using RHS = function<double(double, double)>;

 /**
  * @brief Single Euler step
  * @param f  RHS function f(y, t)
  * @param y  Current value
  * @param t  Current time
  * @param h  Step size
  * @return   y at t+h
  */
 double euler_step(RHS f, double y, double t, double h);

 /**
  * @brief Single Heun (improved Euler / RK2) step
  */
 double heun_step(RHS f, double y, double t, double h);

 /**
  * @brief Single Midpoint method step
  */
 double midpoint_step(RHS f, double y, double t, double h);

 /**
  * @brief Single classical RK4 step
  */
 double rk4_step(RHS f, double y, double t, double h);

 /**
  * @brief Integrate ODE with Euler method over N steps
  * @param f      RHS function
  * @param y0     Initial value
  * @param t_end  End time (start time = 0)
  * @param N      Number of steps
  * @return Final value y(t_end)
  */
 double euler_solve(RHS f, double y0, double t_end, int N);

 /** @brief Integrate ODE with Heun method */
 double heun_solve(RHS f, double y0, double t_end, int N);

 /** @brief Integrate ODE with Midpoint method */
 double midpoint_solve(RHS f, double y0, double t_end, int N);

 /** @brief Integrate ODE with RK4 method */
 double rk4_solve(RHS f, double y0, double t_end, int N);

 /**
  * @brief Compute MSE between numerical solution and exact solution over N+1 points
  * @param f       RHS of ODE
  * @param exact   Exact solution exact(t)
  * @param y0      Initial value
  * @param t_end   End time
  * @param N       Number of steps
  * @param method  0=Euler, 1=Heun, 2=Midpoint, 3=RK4
  * @return MSE
  */
 double ode_mse(RHS f,
     function<double(double)> exact,
     double y0, double t_end, int N,
     int method);