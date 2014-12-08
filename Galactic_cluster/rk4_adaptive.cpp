#include <rk4_adaptive.h>
#include <vec3.h>
#include <cmath>
#include <vector_operations.h>

using std::vector;

RK4_adaptive::RK4_adaptive()
{
}

void RK4_adaptive::RKF_45(Cluster &system, double dt_initial, double total_time)
{
    // Creating vectors
    int n_bodies = system.numberOfBodies();

    vector<vec3> A  = vector<vec3>(2*n_bodies);
    vector<vec3> K1 = vector<vec3>(2*n_bodies);
    vector<vec3> K2 = vector<vec3>(2*n_bodies);
    vector<vec3> K3 = vector<vec3>(2*n_bodies);
    vector<vec3> K4 = vector<vec3>(2*n_bodies);
    vector<vec3> K5 = vector<vec3>(2*n_bodies);
    vector<vec3> K6 = vector<vec3>(2*n_bodies);

    // Setting up A
    for(int i=0;i<n_bodies;i++)
    {
        CelestialBody body = system.bodies[i];
        A[2*i]  = body.position;
        A[2*i+1]= body.velocity;
    }

    double dt = dt_initial;
    double time_elapsed = 0;

    // Maximum allowed relative deviation of RK4 approximation from RK5 approximation
    double tolerance = 1e-6;

    // Minimum and maximum allowed value of timestep
    double mindt = 1e-30;
    double maxdt = 1e-3;

    // Write first step to file for plotting
    system.dumpToFile(time_elapsed);

    // Time loop
    while(time_elapsed<=total_time)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*(time_elapsed)/(total_time));


        // Setting up A
        /* vector A = [position_1, velocity_1, position_2, velocity_2,  ... , position_(n_bodies), velocity_(n_bodies)]
                           (first body)            (second body)                 for all celestial bodies
        */
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            A[2*i]  = body.position;
            A[2*i+1]= body.velocity;
        }

        //----------------------------------------------------------------
        // Calculating K vectors, K1 is equal for RK4 and RKF45
        K1= dAdt(system,A);

        // Returning acceleration values to CelestialBody objects to track acceleration
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            body->acceleration=K1[2*i+1];
        }

        // Calculating K vectors, continued
        K1 = K1*dt;

        //---------------------------------------------------------------------------------------------
        // RK4 approximation using basic RK4 weighting

        K2 = dAdt(system, A+K1*(1/2.0)) * dt;

        K3 = dAdt(system, A+K2*(1/2.0)) * dt;

        K4 = dAdt(system, A+K3        ) * dt;

        vector<vec3> rk4approx = vector<vec3>(2*n_bodies);

        rk4approx = A + ((K1 + K2*2 + K3*2 + K4)*(1/6.0));

        //---------------------------------------------------------------------------------------------
        // RK5 approximation using RKF45 weighting

        K2 = dAdt(system, A + K1*(1/4.0) ) * dt;

        K3 = dAdt(system, A + K1*(3/32.0) + K2*(9/32.0)) * dt;

        K4 = dAdt(system, A + K1*(1932/2197.0) + K2*(-7200/2197.0) + K3*(7296/2197.0) ) * dt;

        K5 = dAdt(system, A + K1*(439/216.0) + K2*(-8) + K3*(3680/513.0) + K4*(-845/4104.0)) * dt;

        K6 = dAdt(system, A + K1*(-8/27.0) + K2*2 + K3*(-3544/2565.0) + K4*(1859/4104.0) + K5*(-11/40.0)) * dt;

        // RK5 approximation
        vector<vec3> rk5approx = vector<vec3>(2*n_bodies);

        rk5approx = A + K1*(16/135.0) + K3*(6656/12825.0) + K4*(28561/56430.0) + K5*(-9/50.0) + K6*(2/55.0);


        //-------------------------------------------------------------------------------------------------------
        // Comparing results of calculations
        double diffMethods = 0;
        double rk5absvalue = 0;
        for (int i=0;i<2*n_bodies;i+=2)
        {
            vec3 diffVector = rk5approx[i]-rk4approx[i];
            diffMethods += fabs(diffVector.x()) + fabs(diffVector.y()) + fabs(diffVector.z());
            rk5absvalue += fabs(rk5approx[i].x()) + fabs(rk5approx[i].y()) + fabs(rk5approx[i].z());
        }
        diffMethods = diffMethods/rk5absvalue;

        if(diffMethods > tolerance && dt>mindt)
        {
            double s = pow(tolerance/(2*diffMethods),1/4.0);
            if(dt*s>mindt)
            {
                dt = dt*s;
                printf("Time step shrunk by a factor %f to %e \n",s,dt);
            }
            else
            {
                dt = dt*0.5;
                printf("Time step halved to %e \n",dt);
            }
        }
        else if(diffMethods < tolerance*0.0001 && dt<maxdt)
        {
            double s = pow(tolerance/(2*diffMethods),1/4.0);
            if(dt*s<=maxdt)
            {
                dt = dt*s;
                printf("Time step increased by a factor %f to %e \n",s,dt);
            }
            else
            {
                dt = dt*2.0;
                printf("Time step doubled to %e \n",dt);
            }
        }
        else
        {
            // Returning new position and velocity values to CelestialBody objects using the accepted RK5 approximation
            for(int i=0;i<n_bodies;i++)
            {
                CelestialBody *body = &system.bodies[i];
                body->position = rk5approx[2*i];
                body->velocity = rk5approx[2*i+1];
            }
            time_elapsed += dt;
            // Write to file for plotting
            system.dumpToFile(time_elapsed);
        }
    }
}

vector<vec3> RK4_adaptive::dAdt(Cluster &system, vector<vec3> A)
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
            double dr_cubed = dr*dr*dr+1e-6;

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
