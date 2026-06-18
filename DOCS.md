# Biblioteka Numeryczna — Dokumentacja

Biblioteka C++14 zawierająca metody numeryczne, bez zewnętrznych zależności. Wszystkie funkcje przyjmują lambdy lub obiekty `std::function`.

---

## Kompilacja

```bash
mkdir build && cd build
cmake ..
make
```

---

## Moduły

- [Równania nieliniowe — znajdowanie pierwiastków](#1-równania-nieliniowe--znajdowanie-pierwiastków)
- [Interpolacja](#2-interpolacja)
- [Solvery RRZ](#3-solvery-rrz)
- [Układy równań liniowych](#4-układy-równań-liniowych)
- [Całkowanie numeryczne](#5-całkowanie-numeryczne)
- [Aproksymacja metodą najmniejszych kwadratów](#6-aproksymacja-metodą-najmniejszych-kwadratów)

---

## 1. Równania nieliniowe — znajdowanie pierwiastków

**Nagłówek:** `include/nonlinear.h`

Szukanie pierwiastków równania f(x) = 0. Wszystkie metody iterują do momentu spełnienia warunku `|b - a| < eps` lub osiągnięcia maksymalnej liczby iteracji `NL_MAX_ITER = 150`. Domyślna tolerancja to `NL_EPS = 1e-9`.

```cpp
#include "nonlinear.h"
auto f  = [](double x){ return x*x - 2.0; };
auto df = [](double x){ return 2.0*x; };
double root = newton_method(f, df, 2.0); // ≈ sqrt(2)
```

---

### `bisection`

```cpp
double bisection(function<double(double)> f, double a, double b, double eps = NL_EPS);
```

Wielokrotnie dzieli przedział `[a, b]` na pół i zachowuje podprzedział, w którym f zmienia znak.

**Wymaga:** f(a) i f(b) mają przeciwne znaki. W przeciwnym razie rzuca `invalid_argument`.

**Zbieżność:** Liniowa — O(h) na iterację, gwarantowana.

| Parametr | Opis |
|----------|------|
| `f` | Ciągła funkcja |
| `a`, `b` | Granice przedziału |
| `eps` | Tolerancja zbieżności |

---

### `regula_falsi`

```cpp
double regula_falsi(function<double(double)> f, double a, double b, double eps = NL_EPS);
```

Podobna do bisekcji, ale dzieli przedział w miejscu przecięcia z osią X prostej przechodzącej przez (a, f(a)) i (b, f(b)), zamiast w połowie.

**Wymaga:** f(a) i f(b) mają przeciwne znaki.

**Zbieżność:** Superliniowa w większości przypadków; może być wolna dla asymetrycznych funkcji.

---

### `newton_method`

```cpp
double newton_method(function<double(double)> f, function<double(double)> df,
                     double x0, double eps = NL_EPS);
```

Iteracja Newtona-Raphsona: `x_{n+1} = x_n - f(x_n) / f'(x_n)`.

**Zbieżność:** Kwadratowa w pobliżu prostego pierwiastka — bardzo szybka.

| Parametr | Opis |
|----------|------|
| `f` | Funkcja |
| `df` | Analityczna pochodna f |
| `x0` | Punkt startowy |
| `eps` | Tolerancja zbieżności |

---

### `newton_numeric`

```cpp
double newton_numeric(function<double(double)> f, double x0,
                      double eps = NL_EPS, double h = 1e-5);
```

Metoda Newtona-Raphsona z pochodną przybliżoną różnicami centralnymi:
`f'(x) ≈ (f(x+h) - f(x-h)) / (2h)`.

Stosować gdy analityczna pochodna jest niedostępna.

| Parametr | Opis |
|----------|------|
| `f` | Funkcja |
| `x0` | Punkt startowy |
| `eps` | Tolerancja zbieżności |
| `h` | Krok różniczkowania numerycznego |

---

### `secant_method`

```cpp
double secant_method(function<double(double)> f, double x0, double x1,
                     double eps = NL_EPS);
```

Przybliża pochodną sieczną przechodzącą przez dwa ostatnie punkty iteracji. Nie wymaga podania pochodnej.

**Zbieżność:** Superliniowa — rząd ≈ 1.618 (złoty podział).

| Parametr | Opis |
|----------|------|
| `f` | Funkcja |
| `x0`, `x1` | Dwa punkty startowe |
| `eps` | Tolerancja zbieżności |

---

### `bisection_iterations`

```cpp
int bisection_iterations(function<double(double)> f, double a, double b, double eps = NL_EPS);
```

Zwraca liczbę iteracji potrzebnych bisekcji do zmniejszenia przedziału poniżej `eps`. Przydatne do porównywania metod.

---

### `newton_iterations`

```cpp
int newton_iterations(function<double(double)> f, function<double(double)> df,
                      double x0, double eps = NL_EPS);
```

Zwraca liczbę iteracji potrzebnych metodzie Newtona-Raphsona do zbieżności z `x0`.

---

## 2. Interpolacja

**Nagłówek:** `include/interpolation.h`

Buduje wielomian przechodzący dokładnie przez zadany zbiór węzłów (xᵢ, yᵢ).

```cpp
#include "interpolation.h"
vector<double> xs = {0.0, 1.0, 2.0};
vector<double> ys = {1.0, 3.0, 7.0};
double val = lagrange_interpolate(1.5, xs, ys);
```

---

### `lagrange_interpolate`

```cpp
double lagrange_interpolate(double x, const vector<double>& xs, const vector<double>& ys);
```

Oblicza wartość wielomianu interpolacyjnego Lagrange'a w punkcie `x`:

```
L(x) = Σ yᵢ · ∏_{j≠i} (x - xⱼ) / (xᵢ - xⱼ)
```

**Złożoność:** O(n²) na obliczenie. Najlepszy dla małej liczby węzłów lub jednorazowych obliczeń.

| Parametr | Opis |
|----------|------|
| `x` | Punkt obliczenia |
| `xs` | Węzły — współrzędne x |
| `ys` | Wartości funkcji w węzłach |

---

### `newton_divided_differences`

```cpp
vector<vector<double>> newton_divided_differences(const vector<double>& xs,
                                                   const vector<double>& ys);
```

Buduje tablicę różnic dzielonych `f[i][j]` potrzebną przez `newton_interpolate`. Wywołać raz, a tablicę wielokrotnie używać dla różnych punktów.

Zwraca górnotrójkątną tablicę 2D.

---

### `newton_interpolate`

```cpp
double newton_interpolate(double x, const vector<double>& xs,
                          const vector<vector<double>>& dd);
```

Oblicza wartość wielomianu Newtona przy użyciu gotowej tablicy różnic dzielonych schematem Hornera:

```
N(x) = f[0][0] + (x-x₀)(f[0][1] + (x-x₁)(f[0][2] + ...))
```

Bardziej numerycznie stabilna niż Lagrange przy wielokrotnych obliczeniach.

| Parametr | Opis |
|----------|------|
| `x` | Punkt obliczenia |
| `xs` | Węzły — współrzędne x |
| `dd` | Tablica z `newton_divided_differences()` |

---

### `interpolation_mse`

```cpp
double interpolation_mse(const vector<double>& xs_all, const vector<double>& ys_all,
                         const vector<double>& xs_nodes, const vector<double>& ys_nodes);
```

Oblicza średni błąd kwadratowy wielomianu interpolacyjnego na pełnym zbiorze danych:

```
MSE = (1/N) · Σ (interpolowana(xᵢ) - yᵢ)²
```

Przydatne do oceny strategii doboru węzłów.

| Parametr | Opis |
|----------|------|
| `xs_all`, `ys_all` | Wszystkie punkty danych |
| `xs_nodes`, `ys_nodes` | Podzbiór użyty jako węzły interpolacji |

---

### `horner_eval`

```cpp
double horner_eval(const vector<double>& a, double x);
```

Oblicza wartość wielomianu w postaci naturalnej `a[0] + a[1]x + a[2]x² + ...` schematem Hornera.

**Złożoność:** O(n) — minimalna liczba mnożeń, numerycznie stabilny.

```cpp
vector<double> a = {1.0, 0.0, 1.0}; // 1 + x²
double v = horner_eval(a, 3.0);      // zwraca 10.0
```

---

## 3. Solvery RRZ

**Nagłówek:** `include/ode.h`

Numeryczne rozwiązywanie skalarnego zagadnienia początkowego:

```
dy/dt = f(y, t),   y(0) = y₀
```

Wszystkie funkcje `_solve` całkują od `t = 0` do `t = t_end` w dokładnie `N` krokach i zwracają wartość końcową.

```cpp
#include "ode.h"
auto f = [](double y, double t){ return -y; }; // dy/dt = -y
double y_end = rk4_solve(f, 1.0, 3.0, 1000);  // ≈ e^(-3)
```

Alias typu `RHS = function<double(double, double)>` jest używany dla prawej strony f(y, t).

---

### `euler_step`

```cpp
double euler_step(RHS f, double y, double t, double h);
```

Wykonuje pojedynczy krok metody Eulera: `y_{n+1} = y_n + h · f(y_n, t_n)`.

**Rząd:** 1 — globalny błąd O(h).

---

### `heun_step`

```cpp
double heun_step(RHS f, double y, double t, double h);
```

Pojedynczy krok metody Heuna (ulepszony Euler / RK2). Prognozuje Eulerem, następnie koryguje:

```
k₁ = f(y, t)
y_nowe = y + h · (k₁ + f(y + h·k₁, t+h)) / 2
```

**Rząd:** 2 — globalny błąd O(h²).

---

### `midpoint_step`

```cpp
double midpoint_step(RHS f, double y, double t, double h);
```

Pojedynczy krok metody punktu środkowego. Oblicza prawą stronę w połowie kroku:

```
mid = y + (h/2) · f(y, t)
y_nowe = y + h · f(mid, t + h/2)
```

**Rząd:** 2 — globalny błąd O(h²).

---

### `rk4_step`

```cpp
double rk4_step(RHS f, double y, double t, double h);
```

Pojedynczy krok klasycznej metody Runge-Kutty 4. rzędu:

```
k₁ = f(y, t)
k₂ = f(y + h/2·k₁, t + h/2)
k₃ = f(y + h/2·k₂, t + h/2)
k₄ = f(y + h·k₃, t + h)
y_nowe = y + (h/6)·(k₁ + 2k₂ + 2k₃ + k₄)
```

**Rząd:** 4 — globalny błąd O(h⁴). Zalecana metoda domyślna.

---

### `euler_solve` / `heun_solve` / `midpoint_solve` / `rk4_solve`

```cpp
double euler_solve   (RHS f, double y0, double t_end, int N);
double heun_solve    (RHS f, double y0, double t_end, int N);
double midpoint_solve(RHS f, double y0, double t_end, int N);
double rk4_solve     (RHS f, double y0, double t_end, int N);
```

Wszystkie cztery mają ten sam interfejs. Całkują `f` od `t=0` do `t=t_end` przy użyciu `N` równomiernych kroków i zwracają `y(t_end)`.

| Parametr | Opis |
|----------|------|
| `f` | Prawa strona f(y, t) |
| `y0` | Wartość początkowa y(0) |
| `t_end` | Czas końcowy |
| `N` | Liczba kroków |

---

### `ode_mse`

```cpp
double ode_mse(RHS f, function<double(double)> exact,
               double y0, double t_end, int N, int method);
```

Oblicza średni błąd kwadratowy między rozwiązaniem numerycznym a dokładnym w `N+1` równomiernie rozmieszczonych punktach czasowych.

| `method` | Użyty solver |
|----------|-------------|
| 0 | Euler |
| 1 | Heun |
| 2 | Punkt środkowy |
| 3 | RK4 |

---

## 4. Układy równań liniowych

**Nagłówek:** `include/linear_systems.h`

Rozwiązuje układy równań liniowych Ax = b.

```cpp
#include "linear_systems.h"
vector<vector<double>> A = {{2, 1}, {5, 3}};
vector<double> b = {8, 19};
auto x = gauss_solve(A, b); // x = {5, -2}
```

---

### `lu_decompose`

```cpp
void lu_decompose(const vector<vector<double>>& A,
                  vector<vector<double>>& L,
                  vector<vector<double>>& U,
                  int n);
```

Rozkłada A = L · U metodą Doolittlea, gdzie L jest jednostkową macierzą dolnotrójkątną (diagonala = 1), a U — górnotrójkątną. Bez wyboru elementu głównego.

**Złożoność:** O(n³).

| Parametr | Opis |
|----------|------|
| `A` | Wejściowa macierz n×n |
| `L` | Wyjście: czynnik dolnotrójkątny |
| `U` | Wyjście: czynnik górnotrójkątny |
| `n` | Rozmiar macierzy |

---

### `forward_substitution`

```cpp
vector<double> forward_substitution(const vector<vector<double>>& L,
                                     const vector<double>& b, int n);
```

Rozwiązuje układ z macierzą dolnotrójkątną Lz = b podstawiając od góry do dołu. Używane wewnętrznie przez `lu_solve`.

**Złożoność:** O(n²).

---

### `backward_substitution`

```cpp
vector<double> backward_substitution(const vector<vector<double>>& U,
                                      const vector<double>& z, int n);
```

Rozwiązuje układ z macierzą górnotrójkątną Ux = z podstawiając od dołu do góry.

Rzuca `runtime_error` gdy U jest osobliwa (element diagonalny bliski zeru).

**Złożoność:** O(n²).

---

### `lu_solve`

```cpp
vector<double> lu_solve(const vector<vector<double>>& A, const vector<double>& b);
```

Rozwiązuje Ax = b w trzech krokach: rozkłada A → L, U; rozwiązuje Lz = b; rozwiązuje Ux = z.

Efektywne gdy ten sam układ rozwiązuje się dla wielu wektorów b (rozkład raz, podstawianie wielokrotnie).

---

### `gauss_solve`

```cpp
vector<double> gauss_solve(vector<vector<double>> A, vector<double> b);
```

Rozwiązuje Ax = b eliminacją Gaussa z częściowym wyborem elementu głównego (zamienia wiersze tak, by na diagonali znajdował się największy element). Bardziej numerycznie stabilna niż `lu_solve` dla źle uwarunkowanych macierzy.

Przyjmuje `A` i `b` przez wartość — oryginały nie są modyfikowane.

**Złożoność:** O(n³).

---

### `analyze_system`

```cpp
string analyze_system(vector<vector<double>> A, vector<double> b, int n);
```

Określa typ układu równań liniowych obliczając rzędy macierzy A i macierzy rozszerzonej [A|b].

| Wartość zwracana | Znaczenie |
|-----------------|-----------|
| `"OZNACZONY"` | Jednoznaczne rozwiązanie — rząd(A) = rząd([A\|b]) = n |
| `"NIEOZNACZONY"` | Nieskończenie wiele rozwiązań — rząd(A) = rząd([A\|b]) < n |
| `"SPRZECZNY"` | Brak rozwiązania — rząd(A) ≠ rząd([A\|b]) |

---

### `residual_norm`

```cpp
double residual_norm(const vector<vector<double>>& A,
                     const vector<double>& x,
                     const vector<double>& b, int n);
```

Oblicza normę euklidesową residuum: `||Ax - b||₂`.

Powinno być bliskie zeru dla dokładnego rozwiązania. Służy do weryfikacji wyników dowolnego solvera.

---

## 5. Całkowanie numeryczne

**Nagłówek:** `include/integration.h`

Przybliża całkę oznaczoną ∫ₐᵇ f(x) dx.

```cpp
#include "integration.h"
auto f = [](double x){ return x * x; };
double I = simpson(f, 0.0, 1.0, 1000); // ≈ 1/3
```

---

### `trapezoidal`

```cpp
double trapezoidal(function<double(double)> f, double a, double b, int m);
```

Złożona reguła trapezów z `m` podprzedziałami o szerokości `h = (b-a)/m`:

```
I ≈ h/2 · [f(a) + 2·f(x₁) + 2·f(x₂) + ... + f(b)]
```

**Błąd:** O(h²) — dokładność drugiego rzędu.

| Parametr | Opis |
|----------|------|
| `f` | Funkcja podcałkowa |
| `a`, `b` | Granice całkowania |
| `m` | Liczba podprzedziałów (≥ 1) |

---

### `simpson`

```cpp
double simpson(function<double(double)> f, double a, double b, int m);
```

Złożona reguła Simpsona. Wymaga parzystego `m` (automatycznie zwiększane o 1 gdy nieparzyste). Używa aproksymacji parabolicznej na każdej parze podprzedziałów:

```
I ≈ h/3 · [f(a) + 4·f(x₁) + 2·f(x₂) + 4·f(x₃) + ... + f(b)]
```

**Błąd:** O(h⁴) — dokładność czwartego rzędu. Dokładna dla wielomianów stopnia ≤ 3.

| Parametr | Opis |
|----------|------|
| `f` | Funkcja podcałkowa |
| `a`, `b` | Granice całkowania |
| `m` | Liczba podprzedziałów (parzysta, ≥ 2) |

---

### `gauss_legendre`

```cpp
double gauss_legendre(function<double(double)> f, double a, double b, int n_nodes);
```

Kwadratura Gaussa-Legendre'a na pojedynczym przedziale `[a, b]` przy użyciu `n_nodes` węzłów. Transformuje całkę na przedział referencyjny [-1, 1] i stosuje preobliczone pary węzeł/waga.

**Dokładna dla wielomianów stopnia ≤ 2·n_nodes − 1.**

| `n_nodes` | Stopień wielomianu całkowanego dokładnie |
|-----------|------------------------------------------|
| 2 | ≤ 3 |
| 3 | ≤ 5 |
| 4 | ≤ 7 |
| 5 | ≤ 9 |

Rzuca `invalid_argument` dla `n_nodes` spoza {2, 3, 4, 5}.

---

### `gauss_legendre_composite`

```cpp
double gauss_legendre_composite(function<double(double)> f,
                                 double a, double b,
                                 int m, int n_nodes);
```

Dzieli `[a, b]` na `m` równych podprzedziałów i stosuje `gauss_legendre` na każdym z nich. Łączy wysoką dokładność Gaussa-Legendre'a z elastycznością metod złożonych dla funkcji silnie zmiennych na przedziale.

| Parametr | Opis |
|----------|------|
| `f` | Funkcja podcałkowa |
| `a`, `b` | Granice całkowania |
| `m` | Liczba podprzedziałów |
| `n_nodes` | Węzły Gaussa na podprzedzial (2–5) |

---

## 6. Aproksymacja metodą najmniejszych kwadratów

**Nagłówek:** `include/approximation.h`

Wyznacza wielomian P(x) = c₀ + c₁x + ... + c_degree·x^degree zadanego stopnia minimalizujący scałkowany błąd kwadratowy ∫ₐᵇ (f(x) − P(x))² dx.

```cpp
#include "approximation.h"
auto f = [](double x){ return (x - 1) * exp(sin(x * x)); };
auto coeffs = least_squares_poly(f, -2.5, 1.0, 10);
double val = eval_approx(coeffs, 0.5);
```

---

### `least_squares_poly`

```cpp
vector<double> least_squares_poly(function<double(double)> f,
                                   double a, double b,
                                   int degree,
                                   int n_int = 1000);
```

Buduje i rozwiązuje równania normalne w celu wyznaczenia współczynników najlepszego wielomianu aproksymującego:

```
A·c = B
A[i][j] = ∫ xⁱ · xʲ dx      (iloczyny skalarne funkcji bazowych)
B[i]    = ∫ f(x) · xⁱ dx    (iloczyny skalarne z aproksymowaną funkcją)
```

Iloczyny skalarne obliczane są regułą Simpsona z `n_int` podprzedziałami. Powstały układ liniowy rozwiązywany jest eliminacją Gaussa.

| Parametr | Opis |
|----------|------|
| `f` | Aproksymowana funkcja |
| `a`, `b` | Przedział aproksymacji |
| `degree` | Stopień wielomianu |
| `n_int` | Podprzedziały całkowania (więcej = dokładniejsze iloczyny skalarne) |

**Zwraca:** Wektor współczynników `c`, gdzie `c[i]` jest współczynnikiem przy xⁱ.

---

### `eval_approx`

```cpp
double eval_approx(const vector<double>& coeffs, double x);
```

Oblicza P(x) = Σ coeffs[i] · xⁱ w zadanym punkcie. Bezpośrednie sumowanie, O(stopień).

---

### `approx_mse`

```cpp
double approx_mse(function<double(double)> f,
                  const vector<double>& coeffs,
                  double a, double b,
                  int steps = 1000);
```

Oblicza scałkowany średni błąd kwadratowy aproksymacji na `[a, b]`:

```
Hn = (1/(b-a)) · ∫ₐᵇ (P(x) - f(x))² dx
```

Przybliżany regułą Simpsona z `steps` podprzedziałami.

---

### `approx_max_error`

```cpp
double approx_max_error(function<double(double)> f,
                        const vector<double>& coeffs,
                        double a, double b,
                        int steps = 1000);
```

Oblicza maksymalny bezwzględny błąd punktowy aproksymacji na `[a, b]`:

```
max |P(xᵢ) - f(xᵢ)|   dla xᵢ równomiernie rozłożonych z krokiem (b-a)/steps
```

---

## Porównanie metod

### Znajdowanie pierwiastków

| Metoda | Wymaga pochodnej | Rząd zbieżności | Uwagi |
|--------|:----------------:|:---------------:|-------|
| Bisekcja | Nie | 1 (liniowy) | Gwarantowana, wolna |
| Regula Falsi | Nie | Superliniowy | Może zwalniać dla asymetrycznych funkcji |
| Newton | Analityczna | 2 (kwadratowy) | Najszybsza przy prostym pierwiastku |
| Newton numeryczny | Nie | ~2 | Gdy pochodna niedostępna |
| Sieczna | Nie | ~1.618 | Dobry wybór ogólny |

### Solvery RRZ

| Metoda | Rząd | Obliczenia prawej strony na krok |
|--------|:----:|:--------------------------------:|
| Euler | 1 | 1 |
| Heun | 2 | 2 |
| Punkt środkowy | 2 | 2 |
| RK4 | 4 | 4 |

### Całkowanie numeryczne

| Metoda | Rząd błędu | Dokładna dla wielomianów do stopnia |
|--------|:----------:|-------------------------------------|
| Trapezy | O(h²) | 1 |
| Simpson | O(h⁴) | 3 |
| Gauss-Legendre (n węzłów) | Bardzo wysoki | 2n − 1 |