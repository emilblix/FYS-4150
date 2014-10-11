#include <iostream>
#include <vec3.h>
#include <cmath>
#include <celestialbody.h>
#include <solarsystem.h>
//#include <runge_kutta_4.h>

using std::cout;
using std::endl;
using std::sqrt;

int main()
{
    // Setting initial solar system and celestial bodies
    const double pi = 3.141592653589793238463;

    SolarSystem solSyst;
    CelestialBody sun(0,0,0,0,1);
    CelestialBody earth(1,0,0,2*pi,3e-6);

    solSyst.addCelestialBody(sun);
    solSyst.addCelestialBody(earth);

    int number_of_years = 2;        // Endpoint of time calculations
    int n_steps = 1000;              // Number of calculation points
    double timestep = (double) number_of_years / (double) n_steps;
    int n_bodies = solSyst.numberOfBodies();

    for(int step=0;step<=n_steps;step++)
    {
//        cout << "step nr "<< step<<endl;
        solSyst.calculateForcesAndEnergy();

        // beregn Velocity og position (for bodies[i=1] til i<numberOfBodies
        for (int i=1;i<n_bodies;i++)
        {
            CelestialBody *body1 = solSyst.bodies[i];
            body1->acceleration = body1->force/body1->mass;
            vec3 da = body1->acceleration*timestep;
            body1->velocity = body1->velocity + da;
            vec3 dv = body1->velocity*timestep;
            body1->position = body1->position + dv;
//        ax = F/m;
//        vx = vx + ax*dt;
//        x = x + v*dt;

        }

//        cout<<"3 : "<< earth.acceleration.x()<<endl;
        solSyst.resetAllForces();
//        cout<<4<< earth.force.x()<<endl;
        solSyst.dumpToFile(timestep, step);
    }


    cout <<  solSyst.numberOfBodies() << endl;
    cout <<  solSyst.bodies.at(1)->mass << endl;
//    solSyst.bodies.at(0)->mass=10;
    cout << sun.mass<<endl;

    return 0;
}


