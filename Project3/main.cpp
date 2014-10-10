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
    SolarSystem solSyst;
    CelestialBody sun(0,0,0,0,1);
    CelestialBody earth(1,0,0,0,3e-6);

    solSyst.addCelestialBody(sun);
    solSyst.addCelestialBody(earth);


    int number_of_years = 2;        // Endpoint of time calculations
    int n_steps = 10;              // Number of calculation points
    double timestep = (double) number_of_years / (double) n_steps;
    int n_bodies = solSyst.numberOfBodies();

    for(int step=0;step<n_steps;step++)
    {
        cout << "step nr "<< step<<endl;
        solSyst.calculateForcesAndEnergy();

        // beregn Velocity og position (for bodies[i=1] til i<numberOfBodies
//        ax = F/m;
//        vx = vx + ax*dt;
//        x = x + v*dt;

        cout<<3<< earth.force.x()<<endl;
        solSyst.resetAllForces();
        cout<<4<< earth.force.x()<<endl;
        solSyst.dumpToFile(timestep, step);
    }


    cout <<  solSyst.numberOfBodies() << endl;
    cout <<  solSyst.bodies.at(1)->mass << endl;
//    solSyst.bodies.at(0)->mass=10;
    cout << sun.mass<<endl;

    return 0;
}


