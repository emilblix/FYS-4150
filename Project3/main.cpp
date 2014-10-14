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
    CelestialBody earth(1/sqrt(2),1/sqrt(2),-2/sqrt(2)*pi,2/sqrt(2)*pi,3e-6);

    solSyst.addCelestialBody(sun);
    solSyst.addCelestialBody(earth);

    float number_of_years = 2;        // Endpoint of time calculations
    int n_steps = 10000;              // Number of calculation points
    double dt = (double) number_of_years / (double) n_steps;
    int n_bodies = solSyst.numberOfBodies();

//    cout<<"timestep = "<<h<<endl;

    for(int step=0;step<=n_steps;step++)
    {
        solSyst.dumpToFile(dt, step);
//        cout << "step nr "<< step<<endl;
        solSyst.calculateForcesAndEnergy();

        // beregn Velocity og position (for bodies[i=1] til i<numberOfBodies
        for (int i=1;i<n_bodies;i++)
        {
            CelestialBody *body1 = solSyst.bodies[i];
            body1->acceleration = body1->force/body1->mass;
            body1->acceleration = body1->acceleration*39.385824;
//            cout<<"body1->position         = "<<body1->position<<endl;
//            cout<<"body1->velocity         = "<<body1->velocity<<endl;
//            cout<<"body1->acceleration     = "<<body1->acceleration<<endl;
            vec3 posnor = body1->position; posnor.normalize();
            vec3 accnor = body1->acceleration; accnor.normalize(); accnor= accnor*(-1);
//            cout<<"body1->position nor     = "<<posnor<<endl;
//            cout<<"body1->acceleration nor = "<<accnor<<endl;
            vec3 nordiff=posnor-accnor;
            double difflength = nordiff.length();
            if(difflength>1e-5)
            {
                cout<<"Difference accnor posnor more than 1e-5"<<endl;
                break;
            }


            vec3 da = body1->acceleration*dt;
            body1->velocity = body1->velocity + da;
            vec3 dv = body1->velocity*dt;
            body1->position = body1->position + dv;
//        ax = F/m;
//        vx = vx + ax*dt;
//        x = x + v*dt;
        }

//        cout<<"earth.force          : "<< earth.force<<endl;
//        cout<<"earth.acceleration a : "<< earth.acceleration<<endl;
//        solSyst.bodies.at(1)->acceleration = solSyst.bodies.at(1)->force/solSyst.bodies.at(1)->mass;
//        cout<<"earth.acceleration b : "<< earth.acceleration<<endl;
        solSyst.resetAllForces();
//        cout<<4<< earth.force.x()<<endl;
    }


//    cout <<  solSyst.numberOfBodies() << endl;
//    cout <<  solSyst.bodies.at(1)->mass << endl;
//    solSyst.bodies.at(0)->mass=10;
//    cout << sun.mass<<endl;

    return 0;
}


