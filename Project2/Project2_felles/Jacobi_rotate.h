#ifndef JACOBI_ROTATE_H
#define JACOBI_ROTATE_H
#include <armadillo>

using namespace arma; // MÅ ha for å kunne bruke mat, vec i funksjoner, fjerner nødvendighet for denne i cpp-fila

void Jacobi_rotate(mat &A, mat &R, int k, int l, int n);
//første ord (int, double, vec, mat) MÅ være samme type som ønsket output (void gir ingen output/return)


#endif // JACOBI_ROTATE_H
