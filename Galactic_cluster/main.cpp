#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>
#include <celestialbody.h>
#include <cluster.h>
#include <rk_4.h>
#include <rk4_adaptive.h>
#include <verlet.h>

using namespace std;


int main()
{
    // INPUT AREA

    // Set name of file to be read. File must contain 7 stats (mass, starting position (x,y,z),
    // starting velocity(vx,vy,vz)) on each line, separated by spaces. Each line counts as one body.
    const char* filename = "../ClusterData/Input/cluster100.txt";

    // Solar system(0) or Galactic cluster(1)?
    int solsystOrCluster = 1;

    // Set integration method. RK4 is 1, adaptive RK4 is 2, Velocity Verlet is 3
    int method = 2;

    // Set endpoint of time calculations and timestep (dt)
    double total_time = 10;
    double timestep = 1e-5;

    //-----------------------------------------------------------------------------------------------
    // Setting initial solar system and celestial bodies
    Cluster astCluster;

    // Add data from file as celestial bodies into the cluster
    ifstream infile(filename);
    double m, x, y, z, vx, vy, vz;
    while (infile >> m >> x >> y >> z >> vx >> vy >> vz)
    {
        CelestialBody body(m, x, y, z, vx, vy, vz);
        astCluster.addCelestialBody(body);
    }

    if(astCluster.numberOfBodies()==0){
        cout << "Error during reading of file. Check filename." << endl;
        return(1);
    }
    cout << "Reading file: " << filename << endl;
    cout << "Number of bodies = "<< astCluster.numberOfBodies()<<endl;


    const double pi = 4*std::atan(1.0); // Pi with double-precision

    //----------------------------------------------------------------------------------

    // FOR SOLAR SYSTEM:
    /* Correction for center of mass, allows input file
       to have the Sun at (0,0,0) with zero velocity and
       removes need to adjust center of mass when adding/removing
       bodies */

    if(solsystOrCluster==0)
    {
        astCluster.gravitationalConstant = 4*pi*pi; // Value of G for solar system calculations

        double xCenter = 0;
        double totMass = 0;
        double totMomentum = 0;
        for(int i=0;i<astCluster.numberOfBodies();i++)
        {
            CelestialBody *body = &astCluster.bodies.at(i);
            xCenter += body->mass*body->position.x();
            totMass += body->mass;
            totMomentum -= body->mass * body->velocity.y();
        }
        xCenter = xCenter/totMass;
        for(int i=0;i<astCluster.numberOfBodies();i++)
        {
            CelestialBody *body = &astCluster.bodies.at(i);
            body->position.set(body->position.x()- xCenter,0,0);
        }
        totMomentum = totMomentum/astCluster.bodies.at(0).mass;
        astCluster.bodies.at(0).velocity.set(0,totMomentum,0);
    }
    else if(solsystOrCluster==1)
    {
        astCluster.gravitationalConstant = 100*pi*pi/(astCluster.numberOfBodies());
    }
    else
    {
        cout << "Error: Invalid value of variable solsystOrCluster" << endl;
        return(1);
    }

    //-----------------------------------------------------------------------------------

    int n_steps = total_time/timestep;              // Number of calculation points

    // Time measurement
    clock_t start, finish;
    start = clock();

    if(method == 1)  // 4th order Runge-Kutta
    {
        cout << "Chosen method is 4th Order Runge-Kutta" << endl;
        RK4::integrateCluster(astCluster,timestep,n_steps);
    }
    else if(method == 2)  // Runge-Kutta-Fehlberg (RKF45)
    {
        cout << "Chosen method is Runge-Kutta-Fehlberg (RKF45)" << endl;
        RK4_adaptive::RKF_45(astCluster,timestep,total_time);
    }
    else if(method == 3)  // Velocity Verlet
    {
        cout << "Chosen method is Velocity Verlet" << endl;
        Verlet::velocityVerlet(astCluster,timestep,n_steps);
    }
    else
    {
        cout << "Error: No chosen integration method" << endl;
        return(1);
    }
    // End timing
    finish = clock();
    double looptime = ((finish-start)/(double) CLOCKS_PER_SEC);

    cout <<endl<< "Simulation complete." << endl;
    cout <<endl<< "Time used in loop: " << looptime << " s" << endl;
    return 0;
}
