#ifndef RK4_H
#define RK4_H
#include <cluster.h>

class RK4
{
public:
    RK4();
    static void integrateCluster(Cluster &system, double dt, int n_steps);
    static std::vector<vec3> dAdt(Cluster &system, std::vector<vec3> A);
};

#endif // RK4_H



