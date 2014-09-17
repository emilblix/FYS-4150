#include <iostream>
#include <cmath>
#include <armadillo>
#include <offdiag.h>
#include "time.h"
using namespace std;
using namespace arma;

double offdiag(mat A, int n)
{


//    vec diagonaltall = linspace(1,n,n);
//    cout<< "diagonaltall= "<<endl<<diagonaltall<<endl;
//    int row, col;
    // Time measurement
    clock_t start, finish;
    start = clock();



    mat X =abs(A);
    X.diag(0).fill(0);
    uword row, col;
//    uword col;
    double max=X.max(row,col);
    finish = clock();
    double time_normal = ((finish-start)/(double) CLOCKS_PER_SEC);
    cout << "Computational time using algorithm: " << time_normal << endl;
/* comp time in secs for values of n:
 n=10^4: 1.34
 n=5*10^3: 0.18
 n=7000: 0.34
 n=15000: 3.11
  */
//    cout << "X= "<<endl<<X<<endl;
    cout<<"max value = "<< max << " , at row "<<row<<" and column "<<col<<endl;


    //    for (int i=0, i<n, i++)
//    {

//    }

    return max;
}
