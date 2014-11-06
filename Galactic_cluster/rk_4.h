#ifndef RK4_H
#define RK4_H
#include <cluster.h>

class RK4
{
public:
    RK4();
//    static void integrateCluster   (Cluster &system, double dt);
//    static std::vector<double> dAdt   (Cluster &system, std::vector<double> A);
//    static std::vector<double> mult(std::vector<double> a, double k);
//    static std::vector<double> add(std::vector<double> a, std::vector<double> b);
    static void integrateClusterVec(Cluster &system, double dt, int n_steps);
    static std::vector<vec3> dAdtVec(Cluster &system, std::vector<vec3> A);
    static std::vector<vec3> multv(std::vector<vec3> a, double k);
    static std::vector<vec3> addv(std::vector<vec3> a, std::vector<vec3> b);
};

#endif // RK4_H



