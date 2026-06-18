# Biblioteka Numeryczna

Biblioteka metod numerycznych napisana w C++14, bez zewnętrznych zależności. Przeznaczona do celów edukacyjnych — pokrywa tematy typowe dla kursu Metod Numerycznych.

## Zawartość

| Moduł | Nagłówek | Opis |
|-------|----------|------|
| Równania nieliniowe | `nonlinear.h` | Bisekcja, Regula Falsi, Newton, Sieczna |
| Interpolacja | `interpolation.h` | Lagrange, Newton (różnice dzielone), Horner |
| Równania różniczkowe | `ode.h` | Euler, Heun, Punkt środkowy, RK4 |
| Układy liniowe | `linear_systems.h` | Dekompozycja LU, Eliminacja Gaussa |
| Całkowanie | `integration.h` | Trapezy, Simpson, Gauss-Legendre |
| Aproksymacja | `approximation.h` | Najmniejsze kwadraty (baza potęgowa) |

## Struktura projektu

```
biblioteka_numeryczna/
├── include/          # Nagłówki publiczne (API biblioteki)
├── src/              # Implementacje
├── tests/            # Testy jednostkowe
├── examples/         # Przykłady użycia
├── DOCS.md           # Szczegółowa dokumentacja wszystkich funkcji
└── CMakeLists.txt    # Konfiguracja budowania
```

## Kompilacja

Wymagania: CMake ≥ 3.10, kompilator C++14 (GCC lub Clang).

```bash
mkdir build && cd build
cmake ..
make
```

Zbudowane pliki:
- Biblioteka statyczna `libnumeryczna.a`
- Testy: `test_nonlinear`, `test_interpolation`, `test_ode`, `test_linear_systems`, `test_integration`, `test_approximation`
- Przykłady: `example_cooling`, `example_roots`

## Uruchomienie testów

```bash
cd build
./test_nonlinear
./test_interpolation
./test_ode
./test_linear_systems
./test_integration
./test_approximation
```

## Szybki start

```cpp
#include "nonlinear.h"
#include "integration.h"
#include "ode.h"

// Znajdź pierwiastek x^2 - 2 = 0
auto f  = [](double x){ return x*x - 2.0; };
auto df = [](double x){ return 2.0*x; };
double root = newton_method(f, df, 2.0);   // ≈ 1.41421356

// Oblicz całkę sin(x) na [0, π]
auto g = [](double x){ return sin(x); };
double I = simpson(g, 0.0, M_PI, 1000);   // ≈ 2.0

// Rozwiąż dy/dt = -y, y(0) = 1
auto rhs = [](double y, double t){ return -y; };
double y_end = rk4_solve(rhs, 1.0, 3.0, 1000);  // ≈ e^(-3)
```

## Dokumentacja

Szczegółowy opis każdej funkcji (parametry, wartości zwracane, złożoność, wzory) znajduje się w pliku [`DOCS.md`](DOCS.md).