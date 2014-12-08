#ifndef VERLET_H
#define VERLET_H
#include <cluster.h>

class Verlet
{
public:
    Verlet();
    static void velocityVerlet(Cluster &system, double dt, int n_steps);
};

#endif // VERLET_H
