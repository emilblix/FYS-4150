#include <iostream>
#include <cmath>
#include <Jacobi_rotate.h>

using namespace std;

void Jacobi_rotate(mat &A, int k, int l, int n)
{
    // Finds values of cos and sin
    double s,c,t;
    double tau = (A(l,l) - A(k,k))/(2*A(k,l));
    if (tau >= 0)
    {
        t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else
    {
        t = -1.0/(-tau + sqrt(1.0 + tau*tau));
    }
    c = 1/sqrt(1 + t * t);
    s = c*t;

    // Saves matrix elements for Jacobi rotation
    double a_ik, a_il; //r_ik, r_il;
    double a_kk = A(k,k);
    double a_ll = A(l,l);

    // Changes the matrix elements with indices k and l
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;

    // Hard-coding the zeros to avoid roundoff errors
    A(k,l) = 0.0; A(l,k) = 0.0;

    // Change remaining elements
    for ( int i = 0; i < n; i++ )
    {
        if ( i != k && i != l )
        {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
    }
    return;
}

