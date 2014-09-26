#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <armadillo>
#include <Jacobi_rotate.h>
#include <offdiag.h>
#include "time.h"

using namespace std;
using namespace arma;


int main()
{

    // Give input values for n, p(max) and omega
    int n, method;
    double pmax, omega;
    cout << "Enter value for n: ";
    cin >> n;
    //    cout << "Enter maximum value for p: ";
    //    cin >> pmax;
//    cout << "Enter value for omega: ";
//    cin >> omega;
    cout << "Enter 1 for single electron, 2 for two electrons without repulsion and 3 for two electrons with repulsion: ";
    cin >> method;

//    n=200;
    pmax=100;
    omega=0.01;

    double h=pmax/n;
    double inverse_hh=1/(h*h);

    // Start matrix A with dimensions (n x n)
    vec rho = linspace(1,n,n)*h;    // rho = (pmin[=0] + i*h), i=0,1,2...N_step
    vec V = vec(n);
    if(method==1){
        V = rho%rho;                                    // For one electron
    }
    else if (method==2){
        V = rho%rho*omega*omega+1/rho;                // For two electrons without repulsion
    }

    mat A=zeros<mat>(n,n);
    A.diag(0) = 2*inverse_hh + V;   // Diagonal = 2/h² + p(i)²
    A.diag(1).fill(-inverse_hh);  // Setting upper and lower tridiagonal as -1
    A.diag(-1).fill(-inverse_hh);



    //        mat R = zeros<mat>(n,n);

    //=======================================================================
    // Eigenvalues for A using Armadillos solver (A must be symmetric)

    // Time measurement
    clock_t start, finish;
    start = clock();

    mat eig_mat_arma = mat(n,n);
    vec eig_val_arma = vec(n);
    eig_sym(eig_val_arma,eig_mat_arma,A);

    // http://stackoverflow.com/questions/24549196/how-to-make-armadillo-work-on-windows

    //    vec eig_val_arma =zeros<vec>(n);
    //    eig_val_arma(0)=3.0;
    //    eig_val_arma(1)=7.0;
    //    eig_val_arma(2)=11.0;

    // End timing for Armadillo calculation
    finish = clock();
    double time_arma = ((finish-start)/(double) CLOCKS_PER_SEC);

    // Saving first eigenvector to file
    mat eigenvec_1 = zeros<mat>(n,2);
    eigenvec_1.col(0)=rho;
    eigenvec_1.col(1)=eig_mat_arma.col(0);

    char *filename = new char[1000];
    sprintf(filename, "Eigenvec_n%04d_o%.3f_p%.1f.dat", n, omega, pmax); cout<<filename<<endl;
    eigenvec_1.save(filename, raw_ascii);


    //=======================================================================
    // Eigenvalues for A using Jacobi rotations

    // Time measurement
    start = clock();

    double tolerance = 10e-5;
    int iterations = 0;
    int maxiter = 1000000;
    double maxnondiag = 1.0;
    int p, q;

    while(maxnondiag > tolerance && iterations <= maxiter)
    {
        maxnondiag = offdiag(A, p, q, n);
        Jacobi_rotate(A, p, q, n);
        iterations++;
    }


    // Eigenvalues must be sorted from smallest to largest
    vec jacobi_eig_val = sort(A.diag(0));



    // End timing for Jacobi method
    finish = clock();
    double time_jacobi = ((finish-start)/(double) CLOCKS_PER_SEC);

    //=======================================================================
    // Printing values

    cout <<endl<< "For p(max) = " << pmax << " and n = " << n << ":" <<endl<<endl;
    cout << "                    Jacobi:     Arma:       Difference: " << endl;
    cout << "1st eigenvalue:     " <<jacobi_eig_val(0)<<"    "<< eig_val_arma(0);
    cout << "    " << jacobi_eig_val(0)-eig_val_arma(0) <<  endl;
    cout << "2nd eigenvalue:     " <<jacobi_eig_val(1)<<"    "<< eig_val_arma(1);
    cout << "    " << jacobi_eig_val(1)-eig_val_arma(1) << endl;
    cout << "3rd eigenvalue:     " <<jacobi_eig_val(2)<<"    "<< eig_val_arma(2);
    cout << "    " << jacobi_eig_val(2)-eig_val_arma(2) << endl;
    cout << endl << "Number of iterations: " << iterations << endl;
    cout << "Time used by Armadillo's solver:     " << time_arma << endl;
    cout << "Time used by Jacobi rotation solver: " << time_jacobi << endl;

    return 0;
}

