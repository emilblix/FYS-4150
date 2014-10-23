#ifndef RK4_H
#define RK4_H
#include <solarsystem.h>

class RK4
{
public:
    RK4();
    void integrateSolarSystem(SolarSystem &system, double dt, int sunCheck);
    std::vector<double> dAdt(SolarSystem &system, std::vector<double> A);
    std::vector<double> mult(std::vector<double> a, double k);
    std::vector<double> add(std::vector<double> a, std::vector<double> b);
};

#endif // RK4_H
