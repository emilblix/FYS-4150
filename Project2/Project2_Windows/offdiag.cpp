#include <iostream>
#include <cmath>
#include <armadillo>
#include <offdiag.h>
#include "time.h"
using namespace std;
using namespace arma;

double offdiag(mat A, int p, int q, int n)
{
    // Time measurement
    clock_t start, finish;
    start = clock();

    double max;


//    mat X =abs(A);
//    X.diag(0).fill(0);
//    uword row, col;
//    double max=X.max(row,col);
// /*comp time in secs for values of n:
// n=5000: 0.18
// n=7000: 0.34
// n=10000: 1.34
// n=15000: 3.11
//  */

//finish = clock();
//double time_normal = ((finish-start)/(double) CLOCKS_PER_SEC);
//cout << "Computational time using .max: " << time_normal << endl;
//cout<<"max value = "<< max << " , at row "<<row<<" and column "<<col<<endl;
//start = clock();
//max = 0;

for (int i=0; i < n; i++)
{
    for (int j= i+1; j<n; j++)
    {
        double aij = fabs(A(i,j));
        if (aij > max)
        {
            max = aij; p = i; q = j;
        }
    }
}

finish = clock();
double time_normal = ((finish-start)/(double) CLOCKS_PER_SEC);
cout << "Computational time using algorithm: " << time_normal << endl;
cout<<"max value = "<< max << " , at row "<<p<<" and column "<<q<<endl;

    return max;
}
