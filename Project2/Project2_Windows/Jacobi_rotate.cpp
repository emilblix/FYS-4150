#include <iostream>
#include <cmath>
#include <armadillo>
#include <Jacobi_rotate.h>
using namespace std;
using namespace arma;

vec Jacobi_rotate(mat A, mat R) // Må ha samme form som "return" (vec for vec, int for int, etc)
{
//    cout <<"R= "<< endl<<R<< endl;
    vec max;
    max= A.col(A.n_cols-2);



    return max;
}

