#include <rk_4.h>
#include <vec3.h>
#include <cmath>
#include <vector_operations.h>

using std::vector;

RK4::RK4()
{
}

void RK4::integrateCluster(Cluster &system, double dt, int n_steps)
{
    // Creating vectors
    int n_bodies = system.numberOfBodies();

    vector<vec3> A  = vector<vec3>(2*n_bodies);
    vector<vec3> K1 = vector<vec3>(2*n_bodies);
    vector<vec3> K2 = vector<vec3>(2*n_bodies);
    vector<vec3> K3 = vector<vec3>(2*n_bodies);
    vector<vec3> K4 = vector<vec3>(2*n_bodies);

    // Time loop
    for(int step=0;step<=n_steps;step++)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

        // Write to file for plotting
        system.dumpToFile(dt*step);

        // Setting up vector A
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            A[2*i]  = body.position;
            A[2*i+1]= body.velocity;
        }

        //----------------------------------------------------------------------------
        // Runge-Kutta 4th order integration, start
        K1= dAdt(system,A);

        // Returning acceleration values to CelestialBody objects to track acceleration
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            body->acceleration=K1[2*i+1];
        }

        // Runge-Kutta 4th order integration, continued

        K1 = K1                         * dt;

        K2 = dAdt(system, A+K1*(1/2.0)) * dt;

        K3 = dAdt(system, A+K2*(1/2.0)) * dt;

        K4 = dAdt(system, A+K3        ) * dt;

        A = A + ((K1 + K2*2.0 + K3*2.0 + K4)*(1/6.0));

        //-----------------------------------------------------------------------------

        // Returning new position and velocity values to CelestialBody objects
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            body->position = A[2*i];
            body->velocity = A[2*i+1];
        }
    }
}

vector<vec3> RK4::dAdt(Cluster &system, vector<vec3> A)
{
    double G = system.gravitationalConstant;
    int n_bodies = system.numberOfBodies();
    vector<vec3> dAdt = vector<vec3>(2*n_bodies);

    // Zeroing the vector dAdt
    for(int i=0;i<2*n_bodies;i++)
    {
        dAdt[i].setToZero();
    }

    // Derivating vector A
    for(int i=0;i<n_bodies;i++)
    {
        // Derivative of position is velocity
        dAdt[2*i]   = A[2*i+1];

        double mass_i = system.bodies.at(i).mass;
        for(int j=i+1; j<n_bodies; j++)
        {
            // Calculating gravitational force between objects
            double mass_j = system.bodies.at(j).mass;
            vec3 dR = A[2*j] - A[2*i];          // dR pointing from body i to body j

            double dr = dR.length();
            double dr_cubed = dr*dr*dr;

            // For body i: a = mass_j/r^3*r
            double forcefactor = mass_j/dr_cubed;
            vec3 valueholder = dR*forcefactor;
            dAdt[2*i+1] = dAdt[2*i+1] + valueholder;

            // For body j: a = mass_i/r^3*r
            forcefactor = mass_i/dr_cubed;
            valueholder = dR*forcefactor;
            dAdt[2*j+1] = dAdt[2*j+1] - valueholder;
        }
        // Multiplying with G at the end
        dAdt[2*i+1] = dAdt[2*i+1]*G;
    }
    return dAdt;
}
