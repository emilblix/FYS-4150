#ifndef VERLET_ADAPTIVE_H
#define VERLET_ADAPTIVE_H
#include <cluster.h>

class Verlet_adaptive
{
public:
    Verlet_adaptive();
    static void velocityVerlet(Cluster &system, double dt_initial, double total_time);
    static void integrateVerlet(Cluster &system, double dt, int n_steps);
    static std::vector<vec3> mult(std::vector<vec3> a, double k);
    static std::vector<vec3> add(std::vector<vec3> a, std::vector<vec3> b);
    static std::vector<vec3> subtract(std::vector<vec3> a, std::vector<vec3> b);
};

#endif // VERLET_ADAPTIVE_H
