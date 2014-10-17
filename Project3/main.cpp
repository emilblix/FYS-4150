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

    // Setting initial solar system and celestial bodies
    const double pi = 4*std::atan(1.0); // Pi with double-presi
    //const double G = 4*pi*pi;

    SolarSystem solSyst;
    CelestialBody sun(0,0,0,0,1,"Sun");
    CelestialBody earth(1,0,0,2*pi,3e-6,"Earth");

    //CelestialBody jupiter(5.2,0,0,2*pi,9.5e-4,"Jupiter");

    solSyst.addCelestialBody(sun);
    solSyst.addCelestialBody(earth);
    //solSyst.addCelestialBody(jupiter);

    RK4 solSystRK;


    float number_of_years = 1;        // Endpoint of time calculations
    double timestep = 1e-4;
    int n_steps = number_of_years/timestep;              // Number of calculation points
    int n_bodies = solSyst.numberOfBodies();


    if(method==0)   // Verlet method
    {
        vector<double> lastAcceleration = vector<double>(2*n_bodies);

        for(int step=0;step<=n_steps;step++)
        {
            printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

            solSyst.dumpToFile(timestep, step);

//            lastAcceleration = avgAcc;
            // previousValue = avgValue

            // Verlet method
            for (int i=1;i<n_bodies;i++) // Inactive
            {
//                CelestialBody *body1 = solSyst.bodies[i];
//                body1->acceleration = body1->force/body1->mass;
//                body1->acceleration = body1->acceleration*G;

//                newValue =
//                avgValue = (newvalue - previousvalue) / 2;
//                body1->position = body1->position + avgvalue*timestep;


//                vec3 da = body1->acceleration*timestep;
//                body1->velocity = body1->velocity + da;
//                vec3 dv = body1->velocity*timestep;
//                body1->position = body1->position + dv;

            }


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
            printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

            solSyst.dumpToFile(timestep, step);

            solSystRK.integrateSolarSystem(solSyst, timestep);

            // Collision test
            int colltest = collisionTest(solSyst);
            if(colltest==1)
            {
                printf(" after %4.1f %% \r", 100.0*((double)step)/((double)n_steps));
                break;
            }

        }
    }
    cout<<endl;
    //    cout <<  solSyst.numberOfBodies() << endl;
    //    cout <<  solSyst.bodies.at(1)->mass << endl;
    //    solSyst.bodies.at(0)->mass=10;
    //    cout << sun.mass<<endl;



    return 0;
}


