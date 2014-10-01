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

    char looptest;
    do{

    // Give input values for n, p(max) and omega
    int n, method;
    double pmax, omega;
    cout << "Enter maximum value for p: ";
    cin >> pmax;
    cout << "Enter value for n_step: ";
    cin >> n;
//    cout << "Enter value for omega: ";
//    cin >> omega;
//    cout << "Enter 1 for single electron, 2 for two electrons without repulsion and 3 for two electrons with repulsion: ";
//    cin >> method;

    omega=0.01;
    method=1;

    double h=pmax/n;
    double inverse_hh=1/(h*h);
    int n_short=n-1;            // Since the matrix is of dimensions n-1 x n-1

    // Start matrix A with dimensions (n-1 x n-1)
    vec rho = linspace(1,n_short,n_short)*h;    // rho = (pmin[=0] + i*h), i=0,1,2...n-1
    vec V = vec(n_short);
    if(method==1){
        V = rho%rho;                                    // For one electron
    }
    else if (method==2){
        V = rho%rho*omega*omega+1/rho;                // For two electrons without repulsion
    }

    mat A=zeros<mat>(n_short,n_short);
    A.diag(0) = 2*inverse_hh + V;   // Diagonal = 2/h² + p(i)²
    A.diag(1).fill(-inverse_hh);  // Setting upper and lower tridiagonal as -1
    A.diag(-1).fill(-inverse_hh);

    // Setting up the eigenvector matrix
    mat R= zeros<mat>(n_short,n_short);
    R.diag(0).fill(1);

    //=======================================================================
    // Eigenvalues for A using Armadillos solver (A must be symmetric)

    // Time measurement
    clock_t start, finish;
    start = clock();

//    mat eig_mat_arma = mat(n_short,n_short);
//    vec eig_val_arma = vec(n_short);
//    eig_sym(eig_val_arma,eig_mat_arma,A);

        vec eig_val_arma =zeros<vec>(n_short);
        eig_val_arma(0)=3.0;
        eig_val_arma(1)=7.0;
        eig_val_arma(2)=11.0;

    // End timing for Armadillo calculation
    finish = clock();
    double time_arma = ((finish-start)/(double) CLOCKS_PER_SEC);

//    // Saving first eigenvector to file
//    mat eigenvec_1 = zeros<mat>(n+1,2);

//    for(int i=1;i<n;i++)
//    {
//        eigenvec_1(i,0)=rho(i-1);
//        eigenvec_1(i,1)=eig_mat_arma(i-1,0);
//    }

//    char *filename = new char[1000];
//    sprintf(filename, "Eigenvec_n%04d_o%.3f_p%.1f.dat", n, omega, pmax);
//    cout<<filename<<endl;
//    eigenvec_1.save(filename, raw_ascii);


    //=======================================================================
    // Eigenvalues for A using Jacobi rotations

    // Time measurement
    start = clock();

    double tolerance = 10e-6;
    int iterations = 0;
    int maxiter = 1000000;
    double maxnondiag = 1.0;
    int p, q;

    while(maxnondiag > tolerance && iterations <= maxiter)
    {
        maxnondiag = offdiag(A, p, q, n_short);
        Jacobi_rotate(A, R, p, q, n_short);
        iterations++;
    }

    // Eigenvalues must be sorted from smallest to largest
    vec jacobi_eig_val = sort(A.diag(0));

    // End timing for Jacobi method
    finish = clock();
    double time_jacobi = ((finish-start)/(double) CLOCKS_PER_SEC);

    //=======================================================================
    // Save first eigenvector to file

        mat eigenvec_1 = zeros<mat>(n+1,2);

        for(int i=1;i<n;i++)
        {
            eigenvec_1(i,0)=rho(i-1);
            eigenvec_1(i,1)=R(i-1,0);
        }

        char *filename = new char[1000];
        sprintf(filename, "Eigenvec_n%04d_o%.3f_p%.1f.dat", n, omega, pmax);
        cout<<filename<<endl;
        eigenvec_1.save(filename, raw_ascii);

    //=======================================================================
    // Printing values

    cout <<endl<< "For p(max) = " << pmax << " and n_step = " << n << ":" <<endl<<endl;
    cout << "                    Jacobi:    Arma:   Difference: " << endl;
    cout << "1st eigenvalue:     " <<jacobi_eig_val(0)<<"    "<< eig_val_arma(0);
    cout << "       " << jacobi_eig_val(0)-eig_val_arma(0) <<  endl;
    cout << "2nd eigenvalue:     " <<jacobi_eig_val(1)<<"    "<< eig_val_arma(1);
    cout << "       " << jacobi_eig_val(1)-eig_val_arma(1) << endl;
    cout << "3rd eigenvalue:     " <<jacobi_eig_val(2)<<"    "<< eig_val_arma(2);
    cout << "      " << jacobi_eig_val(2)-eig_val_arma(2) << endl;
    cout << endl << "Number of iterations: " << iterations << endl;
    cout << "Time used by Armadillo's solver:     " << time_arma << endl;
    cout << "Time used by Jacobi rotation solver: " << time_jacobi << endl;


    //=======================================================================
    // Rerun loop

    cout <<endl<< "New calc? (y/n):";
    cin >> looptest;
    cout << endl;
    }while(looptest!='n');

    return 0;
}

