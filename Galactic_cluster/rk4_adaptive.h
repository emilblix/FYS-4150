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
};

#endif // RK4_ADAPTIVE_H
