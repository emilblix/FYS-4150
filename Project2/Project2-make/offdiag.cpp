#include <iostream>
#include <cmath>
#include <offdiag.h>
#include "time.h"

using namespace std;

double offdiag(const mat &A, int &k, int &l, int n)
{
    double max = 0;

    // Sweeps upper triangle of matrix
    for (int i=0; i < n; i++)
    {
        for (int j= i+1; j<n; j++)
        {
            if (fabs(A(i,j)) > max)
            {
                // Replaces max value if current element is larger
                max = fabs(A(i,j)); k = i; l = j;
            }
        }
    }
    return max;
}
