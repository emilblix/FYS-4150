#ifndef OFFDIAG_H
#define OFFDIAG_H
#include <armadillo>

using namespace arma; // MÅ ha for å kunne bruke mat, vec i funksjoner, fjerner nødvendighet for denne i cpp-fila

double offdiag(const mat &A, int &k, int &l, int n);



#endif // OFFDIAG_H
