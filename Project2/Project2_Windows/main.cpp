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
    cout << "Enter value for p(max): ";
    cin >> pmax;
    cout << "Enter value for n: ";
    cin >> n;

    double h,h_neg2;
    h=pmax/n;
    h_neg2=1/(h*h);

    // Sette opp matriser A og R av dimensjon n x n
    mat A = randn<mat>(n,n);

//        mat A=zeros<mat>(n,n);
//        A.diag(0).fill(2);
//        A.diag(1).fill(-1);  // Setting upper and lower tridiagonal as -1
//        A.diag(-1).fill(-1);
//        A=A*h_neg2;

    mat R = zeros<mat>(n,n);
    vec a1 = A.diag(1);


//    cout << "A= "<<endl<<A<<endl;
//    cout << "A.diag(1)= "<<endl<<a1<<endl;
//    cout << "max A.diag(1)= "<< max(abs(a1))<<endl;
//    cout << "return fra Jacobi_rotate(p,q) er: " << Jacobi_rotate(A,R)<< endl;
    cout<< "offdiagtest: "<<offdiag(A,n)<<endl;



//    double tolerance = 10e-8;
//    int iterations = 0;
//    int maxiter = 1000000;
//    double maxnondiag = 1.0;




//    while(maxnondiag > tolerance && iterations <= maxiter)
//    {
//        int p, q;
//        maxnondiag = offdiag(A, p, q, n);
//        Jacobi_rotate(A, R, p, q, n);
//        iterations++;
//    }



    return 0;
}

