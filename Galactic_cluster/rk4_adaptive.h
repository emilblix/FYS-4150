#ifndef RK4_ADAPTIVE_H
#define RK4_ADAPTIVE_H
#include <cluster.h>


//using std::vector;


class RK4_adaptive
{
public:
//    vector<vec3> halfstep;
//    vector<double> hstepmass;

    RK4_adaptive();
    static void integrateClusterAdaptive(Cluster &system, double dt, int n_steps);
    static std::vector<vec3> dAdt(Cluster &system, std::vector<vec3> A);
    static std::vector<vec3> mult(std::vector<vec3> a, double k);
    static std::vector<vec3> add(std::vector<vec3> a, std::vector<vec3> b);
//    static void halfstep(Cluster system, std::vector<vec3> K4, double dt, int step);
//    static std::vector<vec3> dHalfstepdt(Cluster system, std::vector<vec3> HS, std::vector<vec3> K1);
    static void calcHalfstep(Cluster &system, vector<vec3> &K4, double dt, int step);
};

#endif // RK4_ADAPTIVE_H
