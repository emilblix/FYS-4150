#ifndef RK4_ADAPTIVE_H
#define RK4_ADAPTIVE_H
#include <cluster.h>

class RK4_adaptive
{
public:
    RK4_adaptive();
    static void RKF_45(Cluster &system, double dt_initial, double total_time);
    static std::vector<vec3> dAdt(Cluster &system, std::vector<vec3> A);
    static std::vector<vec3> mult(std::vector<vec3> a, double k);
//    static std::vector<vec3> add(std::vector<vec3> a, std::vector<vec3> b);
//    static std::vector<vec3> add3(std::vector<vec3> a, std::vector<vec3> b, std::vector<vec3> c);
//    static std::vector<vec3> subtract(std::vector<vec3> a, std::vector<vec3> b);
};

#endif // RK4_ADAPTIVE_H
