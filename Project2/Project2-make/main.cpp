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
    // Loop for running the program several times without having to restart it
    char looptest;
    do{

        // Give input values for n, p(max) and omega
        int n, method;
        double pmax, omega;
        cout << "Enter 1 for single electron, 2 for two electrons without repulsion and 3 for two electrons with repulsion: ";
        cin >> method;
        cout << "Rho(max) = ";
        cin >> pmax;
        cout << "n_step   = ";
        cin >> n;
        if(method==2 || method==3)      // Omega is not needed for single electron
        {
            cout << "Omega    = ";
            cin >> omega;
        }

        double h=pmax/n;
        double inverse_hh=1/(h*h);
        int n_short=n-1;            // Must use n_short as the matrix is of dimensions n-1 x n-1

        // Initial matrix A with dimensions (n-1 x n-1)

        vec rho = linspace(1,n_short,n_short)*h;    // rho = (pmin[=0] + i*h), i=0,1,2...n-1
        vec V = vec(n_short);
        if(method==1){
            V = rho%rho;                                    // For one electron
        }
        else if (method==2){
            V = rho%rho*omega*omega;                        // For two electrons without repulsion
        }
        else if (method==3){
            V = rho%rho*omega*omega+1/rho;                  // For two electrons with repulsion
        }

        mat A=zeros<mat>(n_short,n_short);
        A.diag(0) = 2*inverse_hh + V; // Diagonal = 2/h² + p(i)²
        A.diag(1).fill(-inverse_hh);  // Setting upper and lower tridiagonal as -1/h²
        A.diag(-1).fill(-inverse_hh);

        // Setting up the eigenvector matrix
        mat R= zeros<mat>(n_short,n_short);
        R.diag(0).fill(1);

        //========================================================================
        // Eigenvalues for A using Armadillos solver (symmetric A required)

        // Time measurement
        clock_t start, finish;
        start = clock();

        // Calculating eigenvalues and eigenvectors
        mat eig_mat_arma = mat(n_short,n_short);
        vec eig_val_arma = vec(n_short);
        eig_sym(eig_val_arma,eig_mat_arma,A);

        // End timing for Armadillo calculation
        finish = clock();
        double time_arma = ((finish-start)/(double) CLOCKS_PER_SEC);

        // Saving first three eigenvectors to file named with value of omega, pmax and n
        mat eigenvec = zeros<mat>(n+1,4);
        eigenvec(n,0)=pmax;
        for(int i=1;i<n;i++)
        {
            eigenvec(i,0)=rho(i-1);
            eigenvec(i,1)=eig_mat_arma(i-1,0);
            eigenvec(i,2)=eig_mat_arma(i-1,1);
            eigenvec(i,3)=eig_mat_arma(i-1,2);
        }
        char *filename = new char[1000];
        sprintf(filename, "Eigenvec_w%.2f_p%.1f_n%03d.dat", omega, pmax, n);
        eigenvec.save(filename, raw_ascii);


        //========================================================================
        // Eigenvalues for A using Jacobi rotations

        // Time measurement
        start = clock();

        int iterations = 0;
        double tolerance = 10e-6;
        int maxiter = 1000000;
        double maxnondiag = 1.0;
        int p, q;

        while(maxnondiag > tolerance && iterations <= maxiter)
        {
            maxnondiag = offdiag(A, p, q, n_short);
            Jacobi_rotate(A, p, q, n_short);
            iterations++;
        }

        // Eigenvalues must be sorted from smallest to largest
        vec jacobi_eig_val = sort(A.diag(0));

        // End timing for Jacobi method
        finish = clock();
        double time_jacobi = ((finish-start)/(double) CLOCKS_PER_SEC);

        //========================================================================
        // Printing values in a table-like formation

        cout <<endl;
        cout << "                    Jacobi:    Armadillo:   Difference (absolute): " << endl;
        cout << "1st eigenvalue:     " <<jacobi_eig_val(0)<<"    "<< eig_val_arma(0);
        cout << "      " << fabs(jacobi_eig_val(0)-eig_val_arma(0)) <<  endl;
        cout << "2nd eigenvalue:     " <<jacobi_eig_val(1)<<"    "<< eig_val_arma(1);
        cout << "      " << fabs(jacobi_eig_val(1)-eig_val_arma(1)) << endl;
        cout << "3rd eigenvalue:     " <<jacobi_eig_val(2)<<"    "<< eig_val_arma(2);
        cout << "      " << fabs(jacobi_eig_val(2)-eig_val_arma(2)) << endl;
        cout << endl << "Number of iterations: " << iterations << endl;
        cout << "Time used by Armadillo's solver:     " << time_arma << endl;
        cout << "Time used by Jacobi rotation solver: " << time_jacobi << endl;

        //========================================================================
        // Prompt to restart or end program

        cout <<endl<< "New calc? (y/n):";
        cin >> looptest;
        cout <<endl<<"***************************************************************"<<endl<<endl;
    }while(looptest!='n');

    return 0;
}

