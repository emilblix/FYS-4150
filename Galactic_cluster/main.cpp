#include <iostream>
#include <fstream>
#include <celestialbody.h>
#include <cluster.h>
#include <rk_4.h>
#include <verlet.h>
#include <cmath>
#include "time.h"

using namespace std;

int main()
{
    // INPUT AREA

    // Set name of file to be read. File must contain 7 stats (mass, starting position (x,y,z),
    // starting velocity(vx,vy,vz)) on each line, separated by spaces. Each line counts as one body.
    const char* filename = "./Data/solsyst.txt";

    // Set integration method. RK4 is 0, Velocity Verlet is 1, Verlet is 2
    int method = 1;

    // Set endpoint of time calculations and timestep (dt)
    float number_of_years = 50;
    double timestep = 1e-4;

    //-----------------------------------------------------------------------------------------------
    // Setting initial solar system and celestial bodies
    Cluster astCluster;

    // Add data from file as celestial bodies into the cluster
    std::ifstream infile(filename);
    double m, x, y, z, vx, vy, vz;
    while (infile >> m >> x >> y >> z >> vx >> vy >> vz)
    {
        CelestialBody body(m, x, y, z, vx, vy, vz);
        astCluster.addCelestialBody(body);
    }

    cout << "Reading file " << filename << endl;
    cout << "Number of bodies = "<< astCluster.numberOfBodies()<<endl;
    //        for(int i=0;i<astCluster.numberOfBodies();i++)
    //        {
    //            cout << astCluster.bodies.at(i).velocity << endl;
    //        }


    //const double pi = 4*std::atan(1.0); // Pi with double-precision
    //const double G = 4*pi*pi;
    int n_steps = number_of_years/timestep;              // Number of calculation points

    // Time measurement
    clock_t start, finish;
    start = clock();

    cout << "Chosen method is ";
    if(method == 0)  // 4th order Runge-Kutta
    {
        cout << "4th Order Runge-Kutta" << endl;
        RK4::integrateCluster(astCluster,timestep,n_steps);
    }

    else if(method == 1)  // Velocity Verlet
    {
        cout << "Velocity Verlet" << endl;
        Verlet::velocityVerlet(astCluster,timestep,n_steps);
    }

    else if(method == 2)  // Basic Störmer-Verlet
    {
        cout << "Basic Störmer-Verlet" << endl;
        Verlet::integrateVerlet(astCluster,timestep,n_steps);
    }

    // End timing
    finish = clock();
    double looptime = ((finish-start)/(double) CLOCKS_PER_SEC);


    cout <<endl<< "Simulation complete." << endl;
    cout <<endl<< "Time used  in loop: " << looptime << " s" << endl;
    return 0;
}



/* i matlab: dlmwrite */
