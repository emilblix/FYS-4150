#ifndef RK4_H
#define RK4_H
#include <cluster.h>

class RK4
{
public:
    RK4();
    static void integrateCluster(Cluster &system, double dt);
    static std::vector<double> dAdt(Cluster &system, std::vector<double> A);
    static std::vector<double> mult(std::vector<double> a, double k);
    static std::vector<double> add(std::vector<double> a, std::vector<double> b);
};

#endif // RK4_H



