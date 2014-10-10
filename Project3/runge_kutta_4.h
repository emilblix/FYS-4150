#ifndef RUNGE_KUTTA_4_H
#define RUNGE_KUTTA_4_H
#include <armadillo>

using namespace arma; // MÅ ha for å kunne bruke mat, vec i funksjoner, fjerner nødvendighet for denne i cpp-fila

void runge_kutta_4();
//første ord (int, double, vec, mat) MÅ være samme type som ønsket output (void gir ingen output/return)

#endif // RUNGE_KUTTA_4_H
