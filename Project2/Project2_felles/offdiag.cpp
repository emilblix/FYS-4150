#include <iostream>
#include <cmath>
#include <offdiag.h>
#include "time.h"

using namespace std;

double offdiag(const mat &A, int &k, int &l, int n)
{
//    // Time measurement
//    clock_t start, finish;
//    start = clock();

    double max = 0;

    for (int i=0; i < n; i++)
    {
        for (int j= i+1; j<n; j++)
        {
            if (fabs(A(i,j)) > max)
            {
                max = fabs(A(i,j)); k = i; l = j;
            }
        }
    }

//    finish = clock();
//    double time_normal = ((finish-start)/(double) CLOCKS_PER_SEC);
//    cout << "Computational time using algorithm: " << time_normal << endl;
//    cout<<"max value = "<< max << " , at row "<<k<<" and column "<<l<<endl;

    return max;
}
