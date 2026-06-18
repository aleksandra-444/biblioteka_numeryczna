# Testy jednostkowe

Każdy plik testuje jeden moduł biblioteki. Testy nie wymagają zewnętrznych frameworków — korzystają z pomocniczej funkcji `check()` zdefiniowanej lokalnie.

## Pliki

| Plik | Testowany moduł | Liczba funkcji testowych |
|------|-----------------|:------------------------:|
| `test_nonlinear.cpp` | Równania nieliniowe | 6 |
| `test_interpolation.cpp` | Interpolacja | 6 |
| `test_ode.cpp` | Solvery RRZ | 7 |
| `test_linear_systems.cpp` | Układy liniowe | 7 |
| `test_intergration.cpp` | Całkowanie numeryczne | 4 |
| `test_approximation.cpp` | Aproksymacja | 5 |

## Uruchomienie przez CMake (zalecane)

```bash
mkdir -p ../build && cd ../build
cmake ..
make
./test_nonlinear
./test_interpolation
./test_ode
./test_linear_systems
./test_integration
./test_approximation
```

## Kompilacja ręczna (pojedynczy test)

```bash
# Przykład: test równań nieliniowych
g++ -std=c++17 -I../include test_nonlinear.cpp ../src/nonlinear.cpp -o test_nl
./test_nl

# Przykład: test aproksymacji (wymaga dwóch plików źródłowych)
g++ -std=c++17 -I../include test_approximation.cpp \
    ../src/approximation.cpp ../src/linear_systems.cpp -o test_approx
./test_approx
```

## Format wyników

```
======================================
  TESTY: rownienia nieliniowe
======================================

[bisection]:
  [PASS] sqrt(2): blad bisekcji < 1e-9
  [PASS] sqrt(3): blad bisekcji < 1e-9
  [PASS] throws when f(a)*f(b) >= 0
  ...

--------------------------------------
Wyniki: 18 zaliczone, 0 niezaliczone
```

## Co sprawdzają testy

- **Poprawność** — wyniki porównywane z dokładnymi wartościami analitycznymi
- **Dokładność** — błędy mieszczą się w oczekiwanych tolerancjach
- **Obsługa błędów** — nieprawidłowe dane wejściowe rzucają wyjątki
- **Porównanie metod** — weryfikacja kolejności dokładności (np. RK4 > Heun > Euler)
- **Przypadki brzegowe** — wielokrotne pierwiastki, macierze diagonalne, układy sprzeczne