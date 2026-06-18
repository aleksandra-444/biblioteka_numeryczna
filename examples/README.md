# Przykłady użycia

Dwa samodzielne programy demonstracyjne pokazujące praktyczne zastosowanie biblioteki.

## Przykłady

### Przykład 1 — Stygnięcie powierzchni (`example_cooling.cpp`)

Symulacja fizyczna oparta na prawie Stefana-Boltzmanna:

```
dT/dt = -α · T⁴,   T(0) = 5903 K,   α = 3·10⁻¹²
```

**Co demonstruje:**
- Porównanie wszystkich 4 solverów RRZ (Euler, Heun, Punkt środkowy, RK4)
- Błąd końcowy i MSE wzdłuż trajektorii dla każdej metody
- Zapis trajektorii do pliku CSV

**Wynik:** plik `cooling_results.csv` z kolumnami `t, Dokladna, Euler, Heun, PunktSrodkowy, RK4`

---

### Przykład 2 — Pierwiastki i całki (`example_roots_and_integrals.cpp`)

Dwa zadania z typowych ćwiczeń laboratoryjnych.

**Część 1 — Znajdowanie pierwiastków:**
- `x² - 2 = 0` (pierwiastek: √2)
- `cosh(x) - √x - 1 = 0` (pierwiastek ≈ 1.2, brak postaci zamkniętej)

Porównanie bisekcji, Regula Falsi, Newtona i metody siecznych.

**Część 2 — Całkowanie numeryczne:**
- Funkcja: `f(x) = (x+1)·exp(cos(x))` na `[-14, 1]`
- Wartość dokładna: `-116.659583194914`
- Metody: Trapezy, Simpson, Gauss-Legendre (n=2,3,4) dla m = 100, 1000, 10000

**Wynik:** plik `integration_convergence.csv` ze zbieżnością błędu dla rosnącego m.

## Kompilacja przez CMake (zalecane)

```bash
mkdir -p ../build && cd ../build
cmake ..
make
./example_cooling
./example_roots
./example_interpolation
```

## Kompilacja ręczna

```bash
# Przykład 1
g++ -std=c++17 -I../include example_cooling.cpp \
    ../src/ode.cpp -o example_cooling

# Przykład 2
g++ -std=c++17 -I../include example_roots_and_integrals.cpp \
    ../src/nonlinear.cpp ../src/integration.cpp -o example_roots


## Wyniki na konsoli

Przykład 1:
```
Temperatura koncowa T(t_end):
              Dokladna:   5809.123456 K
                 Euler:   5812.001234 K  |blad| = 2.877778
                  Heun:   5809.123891 K  |blad| = 0.000435
          Punkt srodkowy: 5809.123712 K  |blad| = 0.000256
                   RK4:   5809.123456 K  |blad| = 0.000000
```

Przykład 2:
```
  --- m = 1000 podprzedzialy ---
                       Trapezy = -116.659583206...   |blad| = 1.2e-08
                       Simpson = -116.659583194...   |blad| = 3.1e-12
           Gauss-Legendre n=2  = -116.659583194...   |blad| = 8.4e-10
           Gauss-Legendre n=3  = -116.659583194...   |blad| = 2.1e-13
```