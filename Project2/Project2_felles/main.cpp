#include <iostream>
#include <cmath>
#include <armadillo>
#include <Jacobi_rotate.h>
#include <offdiag.h>

using namespace std;
using namespace arma;

int main()
{
    // Gi inputverdi for n og p(max)
    int n;
    double pmax;
    //    cout << "Enter maximum value for p: ";
    //    cin >> pmax;
    //    cout << "Enter value for n: ";
    //    cin >> n;
    n=200; pmax=100;

    double h=pmax/n;
    double inverse_hh=1/(h*h);

    // Sette opp startmatrise A av dimensjon n x n
    vec rho = linspace(1,n,n)*h;    // rho = (pmin[=0] + i*h), i=0,1,2...N_step
    vec V   = rho%rho;              // Element-wise multiplication

    mat A=zeros<mat>(n,n);
    A.diag(0) = 2*inverse_hh + V;   // Diagonal = 2/h² + p(i)²
    A.diag(1).fill(-inverse_hh);  // Setting upper and lower tridiagonal as -1
    A.diag(-1).fill(-inverse_hh);


    //        mat R = zeros<mat>(n,n);


    // Eigenvalues for A using Armadillos solver (A must be symmetric)
    vec eig_val_arma = eig_sym(A);

    double tolerance = 10e-5;
    int iterations = 0;
    int maxiter = 1000000;
    double maxnondiag = 1.0;

    while(maxnondiag > tolerance && iterations <= maxiter)
    {
        int p, q;
        maxnondiag = offdiag(A, p, q, n);
        Jacobi_rotate(A, p, q, n);
        iterations++;
    }


    // Eigenvalues must be sorted from smallest to largest
    vec jacobi_eig_val = sort(A.diag(0));

    cout << "For p(max) = " << pmax << " and n = " << n << ":" <<endl<<endl;
    cout << "                    Jacobi:    " << "Arma: " << endl;
    cout << "1st eigenvalue:     " <<jacobi_eig_val(0)<<"    "<< eig_val_arma(0) << endl;
    cout << "2nd eigenvalue:     " <<jacobi_eig_val(1)<<"    "<< eig_val_arma(1) << endl;
    cout << "3rd eigenvalue:     " <<jacobi_eig_val(2)<<"    "<< eig_val_arma(2) << endl;
    cout << "4th eigenvalue:     " <<jacobi_eig_val(3)<<"    "<< eig_val_arma(3) << endl;
    cout << "5th eigenvalue:     " <<jacobi_eig_val(4)<<"    "<< eig_val_arma(4) << endl;
    cout << "6th eigenvalue:     " <<jacobi_eig_val(5)<<"    "<< eig_val_arma(5) << endl;
    cout << endl << "Difference 1st value: = " << jacobi_eig_val(0)-eig_val_arma(0) << ", diff 2nd value= " << jacobi_eig_val(1)-eig_val_arma(1) << endl;
    cout << "Number of runs through the program loop: " << iterations << endl;


    return 0;
}

