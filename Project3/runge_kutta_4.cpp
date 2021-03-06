#include <runge_kutta_4.h>
#include <vec3.h>
#include <cmath>

using std::vector;

RK4::RK4()
{
}

void RK4::integrateSolarSystem(SolarSystem &system, double dt, int sunCheck)
{
    // Creating vectors
    int n_bodies = system.numberOfBodies();

    std::vector<double> A  = std::vector<double>(4*n_bodies);
    std::vector<double> K1 = std::vector<double>(4*n_bodies);
    std::vector<double> K2 = std::vector<double>(4*n_bodies);
    std::vector<double> K3 = std::vector<double>(4*n_bodies);
    std::vector<double> K4 = std::vector<double>(4*n_bodies);

    // Setting up A
    /* vector A = [x1, y1, vx1, vy1,    (first body)
                   x2, y2, vx2, vy2,    (second body)
                   ... for all celestial bodies
    */
    for(int i=0;i<n_bodies;i++)
    {
        CelestialBody *body = system.bodies[i];
        A[4*i]  = body->position.x();
        A[4*i+1]= body->position.y();
        A[4*i+2]= body->velocity.x();
        A[4*i+3]= body->velocity.y();
    }

    // Runge-Kutta 4th order integration, start
    K1= dAdt(system,A);

    // Returning acceleration values to CelestialBody objects
    for(int i=0;i<n_bodies;i++)
    {
        CelestialBody *body = system.bodies[i];
        body->acceleration.set(K1[4*i+2], K1[4*i+3], 0);
    }

    // Runge-Kutta 4th order integration, continued
    K1 = mult(K1,dt);
    K2 = dAdt(system,add(A,mult(K1,1/2.0))); K2 = mult(K2,dt);
    K3 = dAdt(system,add(A,mult(K2,1/2.0))); K3 = mult(K3,dt);
    K4 = dAdt(system,add(A,K3))            ; K4 = mult(K4,dt);

    // Combining (K1 + 2*K2 + 2*K3 + K4)/6 into K1
    K1 = add(K1,K4);
    K1 = add(K1,mult(K2,2));
    K1 = add(K1,mult(K3,2));
    K1 = mult(K1,1/6.0);
    A = add(A,K1);

    // Returning new position and velocity values to CelestialBody objects
    for(int i=sunCheck;i<n_bodies;i++)  // For stationary Sun, set startpoint i=1
    {
        CelestialBody *body = system.bodies[i];
        body->position.set(A[4*i],   A[4*i+1], 0);
        body->velocity.set(A[4*i+2], A[4*i+3], 0);
    }
}


std::vector<double> RK4::dAdt(SolarSystem &system, std::vector<double> A)
{
    double pi = 4*std::atan(1.0);
    double G = 4*pi*pi;
    int n_bodies = system.numberOfBodies();
    std::vector<double> dAdt = std::vector<double>(4*n_bodies);

    // Zeroing the vector dAdt
    for(int i=0;i<4*n_bodies;i++)
    {
        dAdt[i]=0;
    }

    // Derivating vector A
    for(int i=0;i<n_bodies;i++)
    {
        // Derivative of position is velocity
        dAdt[4*i] = A[4*i+2];
        dAdt[4*i+1] = A[4*i+3];

        double mass_i = system.bodies.at(i)->mass;
        for(int j=i+1; j<n_bodies; j++)
        {
            // Calculating gravitational force between objects
            double mass_j = system.bodies.at(j)->mass;
            double dX = A[4*j]-A[4*i];    // dX and dY "pointing" from body1 to body2
            double dY = A[4*j+1]-A[4*i+1];
            double dr = std::sqrt(dX*dX+dY*dY);

            double forcefactor = G*mass_i*mass_j/(dr*dr*dr);

            dAdt[4*i+2] += dX*forcefactor/mass_i;
            dAdt[4*i+3] += dY*forcefactor/mass_i;

            dAdt[4*j+2] -= dX*forcefactor/mass_j;
            dAdt[4*j+3] -= dY*forcefactor/mass_j;
        }
    }
return dAdt;
}


std::vector<double> RK4::mult(std::vector<double> a, double k) {
    for (unsigned int i=0; i < a.size(); i++) {
        a[i] = a[i]*k;
    }
    return a;
}

std::vector<double> RK4::add(std::vector<double> a, std::vector<double> b) {

    if (a.size() != b.size()) {
        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
        return a;
    }
    std::vector<double> c = std::vector<double>(a.size());
    for (unsigned int i=0; i < a.size(); i++) {
        c[i] = a[i] + b[i];
    }
    return c;
}
