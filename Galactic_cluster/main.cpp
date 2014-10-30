#include <iostream>
#include <fstream>
#include <celestialbody.h>
#include <cluster.h>
#include <rk_4.h>
#include <cmath>

std::ifstream infile("S_E_M.txt");
using namespace std;



int main()
{


    // Set endpoint of time calculations and timestep (dt)
    float number_of_years = 1;
    double timestep = 1e-4;

    // Setting initial solar system and celestial bodies
    const double pi = 4*std::atan(1.0); // Pi with double-precision
//    const double G = 4*pi*pi;

    Cluster astCluster;

    double m, x, y, z, vx, vy, vz;
    while (infile >> m >> x >> y >> z >> vx >> vy >> vz)
    {
        CelestialBody body(m, x, y, z, vx, vy, vz);
        astCluster.addCelestialBody(body);
    }

    cout << "Number of bodies = "<< astCluster.numberOfBodies()<<endl;
    for(int i=0;i<astCluster.numberOfBodies();i++)
    {
        cout << astCluster.bodies.at(i).velocity<<endl;
    }

    int n_steps = number_of_years/timestep;              // Number of calculation points
//    int n_bodies = astCluster.numberOfBodies();

    for(int step=0;step<=n_steps;step++)        // RK4
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

        // Write to file for plotting
        astCluster.dumpToFile(timestep, step);

        RK4::integrateCluster(astCluster,timestep);




    }

    cout <<endl<< "Success!" << endl;
    return 0;
}



/* i matlab: dlmwrite */
