#pragma once
/**
 * @file linear_systems.h
 * @brief Rozwiązywanie układów równań liniowych: dekompozycja LU, eliminacja Gaussa
 *
 * Udostępnia dekompozycję LU (Doolittle), podstawianie wprzód/wstecz
 * oraz eliminację Gaussa z częściowym wyborem elementu głównego.
 *
 * Przykład:
 *   vector<vector<double>> A = {{2,1},{5,3}};
 *   vector<double> b = {8, 19};
 *   auto x = gauss_solve(A, b);
 *   // x = {5, -2}
 */

#include <vector>
#include <string>
using std::vector;
using std::string;

/**
 * @brief Dekompozycja LU (metoda Doolittlea, bez wyboru elementu głównego)
 * @param A Wejściowa macierz kwadratowa (n x n)
 * @param L Wyjściowa macierz dolnotrójkątna (L[i][i] = 1)
 * @param U Wyjściowa macierz górnotrójkątna
 * @param n Rozmiar macierzy
 */
void lu_decompose(const vector<vector<double>>& A,
    vector<vector<double>>& L,
    vector<vector<double>>& U,
    int n);

/**
 * @brief Podstawianie wprzód: rozwiązuje Lz = b
 * @param L Macierz dolnotrójkątna
 * @param b Wektor prawej strony
 * @param n Rozmiar układu
 * @return Wektor rozwiązania z
 */
vector<double> forward_substitution(const vector<vector<double>>& L,
    const vector<double>& b,
    int n);

/**
 * @brief Podstawianie wstecz: rozwiązuje Ux = z
 * @param U Macierz górnotrójkątna
 * @param z Wektor prawej strony
 * @param n Rozmiar układu
 * @return Wektor rozwiązania x
 */
vector<double> backward_substitution(const vector<vector<double>>& U,
    const vector<double>& z,
    int n);

/**
 * @brief Rozwiązuje Ax = b za pomocą dekompozycji LU
 * @param A Macierz kwadratowa (n x n)
 * @param b Wektor prawej strony
 * @return Wektor rozwiązania x
 */
vector<double> lu_solve(const vector<vector<double>>& A,
    const vector<double>& b);

/**
 * @brief Eliminacja Gaussa z częściowym wyborem elementu głównego
 * @param A Macierz kwadratowa (n x n), NIE jest modyfikowana (używana jest kopia)
 * @param b Wektor prawej strony
 * @return Wektor rozwiązania x
 */
vector<double> gauss_solve(vector<vector<double>> A, vector<double> b);

/**
 * @brief Analizuje układ równań liniowych pod kątem spójności i jednoznaczności
 * @param A Macierz współczynników
 * @param b Wektor prawej strony
 * @param n Rozmiar układu
 * @return "OZNACZONY", "NIEOZNACZONY" lub "SPRZECZNY"
 */
string analyze_system(vector<vector<double>> A, vector<double> b, int n);

/**
 * @brief Oblicza normę residuum ||Ax - b||_2
 * @param A Macierz współczynników
 * @param x Wektor rozwiązania
 * @param b Wektor prawej strony
 * @param n Rozmiar układu
 * @return Norma euklidesowa residuum
 */
double residual_norm(const vector<vector<double>>& A,
    const vector<double>& x,
    const vector<double>& b,
    int n);