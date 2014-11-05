#include <iostream>
#include <fstream>
#include <celestialbody.h>
#include <cluster.h>
#include <rk_4.h>
#include <verlet.h>
#include <cmath>
#include "time.h"

using namespace std;

// Set name of file to be read. File must contain 7 stats (mass, starting position,
// starting velocity) on each line, separated by spaces. Each line counts as one body.
std::ifstream infile("S_E_M.txt");

int main()
{
    // Initial values and calculations

    // Set integration method. RK4 is 0, Velocity Verlet is 1, Verlet is 2
    int method = 0;

    // Set endpoint of time calculations and timestep (dt)
    float number_of_years = 5;
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
    //    for(int i=0;i<astCluster.numberOfBodies();i++)
    //    {
    //        cout << astCluster.bodies.at(i).velocity << endl;
    //    }

    //    astCluster.calculateForces();
    //    astCluster.calculateAcceleration();
    //    cout<<endl<<astCluster.bodies.at(1).force<<endl;

    int n_steps = number_of_years/timestep;              // Number of calculation points
    //    int n_bodies = astCluster.numberOfBodies();

    // Time measurement
    clock_t start, finish;
    start = clock();

    if(method == 0)  // 4th order Runge-Kutta
    {
        for(int step=0;step<=n_steps;step++)
        {
            // Print progress bar
            printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

            // Write to file for plotting
            astCluster.dumpToFile(timestep, step);

            RK4::integrateClusterVec(astCluster,timestep);
        }
        // Print for last step
        astCluster.dumpToFile(number_of_years,1);
    }

    if(method == 1)  // Velocity Verlet
    {
        Verlet::velocityVerlet(astCluster,timestep,n_steps);
    }

    if(method == 2)  // Basic Störmer-Verlet
    {
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
