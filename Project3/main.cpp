#include <iostream>
#include <cstdlib>
#include <vec3.h>
#include <cmath>
#include <celestialbody.h>
#include <solarsystem.h>
#include <runge_kutta_4.h>
#include <collisiontest.h>

using std::cout;
using std::endl;
using std::sqrt;
using std::vector;

int main()
{
    // SET METHOD: 0 for Verlet and 1 for Runge-Kutta
    int method = 1;

    // SET SUN MOVEMENT ABILITY: 0 for moving sun, 1 for stationary sun
    int stationarySun = 0;

    // Set endpoint of time calculations and timestep (dt)
    float number_of_years = 400;
    double timestep = 1e-2;

    // Setting initial solar system and celestial bodies
    const double pi = 4*std::atan(1.0); // Pi with double-precision
    const double G = 4*pi*pi;

    // CelestialBody name(x,     y,vx, vy,               mass,    radius,  name);
    CelestialBody sun    (0,     0, 0, 0,                1,       4.64e-3, "Sun");
    CelestialBody mercury(0.39,  0, 0, 2*pi/sqrt(0.39),  1.65e-7, 1.63e-5, "Mercury");
    CelestialBody venus  (0.72,  0, 0, 2*pi/sqrt(0.72),  2.45e-6, 4.03e-5, "Venus");
    CelestialBody earth  (1,     0, 0, 2*pi,             3e-6,    4.25e-5, "Earth");
    CelestialBody mars   (1.52,  0, 0, 2*pi/sqrt(1.52),  3.2e-7,  2.26e-5, "Mars");
    CelestialBody jupiter(5.2,   0, 0, 2*pi/sqrt(5.2),   9.5e-4,  4.66e-4, "Jupiter");
    CelestialBody saturn (9.54,  0, 0, 2*pi/sqrt(9.54),  2.85e-4, 3.88e-4, "Saturn");
    CelestialBody uranus (19.19, 0, 0, 2*pi/sqrt(19.19), 4.35e-5, 1.69e-4, "Uranus");
    CelestialBody neptune(30.07, 0, 0, 2*pi/sqrt(30.07), 5.1e-5,  1.64e-4, "Neptune");
    CelestialBody pluto  (39.48, 0, 0, 2*pi/sqrt(39.48), 6.55e-9, 7.67e-6, "Pluto");

    SolarSystem solSyst;
    solSyst.addCelestialBody(sun);
    solSyst.addCelestialBody(mercury);
    solSyst.addCelestialBody(venus);
    solSyst.addCelestialBody(earth);
    solSyst.addCelestialBody(mars);
    solSyst.addCelestialBody(jupiter);
    solSyst.addCelestialBody(saturn);
    solSyst.addCelestialBody(uranus);
    solSyst.addCelestialBody(neptune);
    solSyst.addCelestialBody(pluto);



    RK4 solSystRK;


    int n_steps = number_of_years/timestep;              // Number of calculation points
    int n_bodies = solSyst.numberOfBodies();


    if(method==0)   // Verlet method
    {
        // Defining vectors for saving acceleration values
        vector<double> newAcceleration = vector<double>(2*n_bodies);
        vector<double> avgAcceleration = vector<double>(2*n_bodies);

        // avgAcceleration must be zero at start
        for(int i=0;i<2*n_bodies;i++)
        {
            avgAcceleration[i]=0;
        }

        for(int step=0;step<=n_steps;step++)
        {
            printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));
            solSyst.dumpToFile(timestep, step);

            // Calculating forces for bodies in solar system
            solSyst.calculateForces();

            for (int i=stationarySun;i<n_bodies;i++)
            {
                CelestialBody *body = solSyst.bodies[i];

                // Fetching forces acting on body
                newAcceleration[2*i]   = G*body->force.x()/body->mass;
                newAcceleration[2*i+1] = G*body->force.y()/body->mass;

                // Setting acceleration as mean value of current step and previous step
                avgAcceleration[2*i]   = (avgAcceleration[2*i]   + newAcceleration[2*i])  /2.0;
                avgAcceleration[2*i+1] = (avgAcceleration[2*i+1] + newAcceleration[2*i+1])/2.0;

                // Returning acceleration to class parameter
                body->acceleration.set(avgAcceleration[2*i], avgAcceleration[2*i+1], 0);

                // Calculating new velocity and returning value to class parameter
                double velocity_x = body->velocity.x()+ avgAcceleration[2*i]  * timestep;
                double velocity_y = body->velocity.y()+ avgAcceleration[2*i+1]* timestep;
                body->velocity.set(velocity_x, velocity_y, 0);

                // Position = previous position + v*dt + 1/2 a*dt^2, returned to calss parameter
                double position_x = body->position.x()+ velocity_x * timestep + ( 0.5 * avgAcceleration[2*i]   * timestep * timestep);
                double position_y = body->position.y()+ velocity_y * timestep + ( 0.5 * avgAcceleration[2*i+1] * timestep * timestep);
                body->position.set(position_x, position_y, 0);

            }

//            cout << "Sum kinetic and potential energy: " << solSyst.totalEnergy()<<endl;

            // Collision test
            int colltest = collisionTest(solSyst);
            if(colltest==1)
            {
                printf(" after %4.1f %% \r", 100.0*((double)step)/((double)n_steps));
                break;
            }
        }
    }

    else if(method==1) // 4th order Runge-Kutta
    {

        for(int step=0;step<=n_steps;step++)
        {
            // Print progress bar
            printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

            // Write to file for plotting
            solSyst.dumpToFile(timestep, step);

            solSystRK.integrateSolarSystem(solSyst, timestep, stationarySun);

//            cout << "Sum kinetic and potential energy: " << solSyst.totalEnergy()<<endl;

            // Collision test
            int colltest = collisionTest(solSyst);
            if(colltest==1)
            {
                printf(" after %4.1f %% \r", 100.0*((double)step)/((double)n_steps));
                break;
            }

        }
    }
    // Write endpoint to file
    solSyst.dumpToFile(timestep,n_steps);
    cout<<endl;


//        cout <<  solSyst.numberOfBodies() << endl;
    //    cout <<  solSyst.bodies.at(1)->mass << endl;
    //    solSyst.bodies.at(0)->mass=10;
    //    cout << sun.mass<<endl;



    return 0;
}


