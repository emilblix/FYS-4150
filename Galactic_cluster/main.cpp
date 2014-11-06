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
    // Set name of file to be read. File must contain 7 stats (mass, starting position,
    // starting velocity) on each line, separated by spaces. Each line counts as one body.
    std::ifstream infile("solsyst.txt");

    // Initial values and calculations

    // Set integration method. RK4 is 0, Velocity Verlet is 1, Verlet is 2
    int method = 0;

    // Set endpoint of time calculations and timestep (dt)
    float number_of_years = 50;
    double timestep = 1e-3;
    //const double pi = 4*std::atan(1.0); // Pi with double-precision
    //const double G = 4*pi*pi;

    // Setting initial solar system and celestial bodies
    Cluster astCluster;

    // Add data from file as celestial bodies into the cluster
    double m, x, y, z, vx, vy, vz;
    while (infile >> m >> x >> y >> z >> vx >> vy >> vz)
    {
        CelestialBody body(m, x, y, z, vx, vy, vz);
        astCluster.addCelestialBody(body);
    }

    cout << "Number of bodies = "<< astCluster.numberOfBodies()<<endl;
    //        for(int i=0;i<astCluster.numberOfBodies();i++)
    //        {
    //            cout << astCluster.bodies.at(i).velocity << endl;
    //        }


    int n_steps = number_of_years/timestep;              // Number of calculation points

    // Time measurement
    clock_t start, finish;
    start = clock();

    if(method == 0)  // 4th order Runge-Kutta
    {
        cout<<"4th Order Runge-Kutta"<<endl;
        RK4::integrateClusterVec(astCluster,timestep,n_steps);
    }

    if(method == 1)  // Velocity Verlet
    {
        cout<<"Velocity Verlet"<<endl;
        Verlet::velocityVerlet(astCluster,timestep,n_steps);
    }

    if(method == 2)  // Basic Störmer-Verlet
    {
        cout<<"Basic Störmer-Verlet"<<endl;
        Verlet::integrateVerlet(astCluster,timestep,n_steps);
    }



    // End timing
    finish = clock();
    double time = ((finish-start)/(double) CLOCKS_PER_SEC);


    cout <<endl<< "Success!" << endl;
    cout <<endl<< "Time used  in loop: " << time << " s" << endl;
    return 0;
}



/* i matlab: dlmwrite */
