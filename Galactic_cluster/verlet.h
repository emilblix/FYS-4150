#ifndef VERLET_H
#define VERLET_H
#include <cluster.h>

class Verlet
{
public:
    Verlet();
    static void velocityVerlet(Cluster &system, double dt, int n_steps);
    static void integrateVerlet(Cluster &system, double dt, int n_steps);
    static std::vector<vec3> mult(std::vector<vec3> a, double k);
    static std::vector<vec3> add(std::vector<vec3> a, std::vector<vec3> b);
};

#endif // VERLET_H
