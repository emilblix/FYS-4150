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

// WARNING!!!!!!!!
// BEFORE DELIVERY:
// SEARCH ENTIRE PROJECT FOR balle
// REMOVE AS NEEDED

/*
//vector<vec3> operator*(vector<vec3> a, double k) {
//    for (unsigned int i=0; i < a.size(); i++)
//    {
//        a[i] = a[i]*k;
//    }
//    return a;
//}

//// Function for adding two std::vector<vec3>, with test for equal length
//vector<vec3> operator+(vector<vec3> a, vector<vec3> b){
//    if (a.size() != b.size())
//    {
//        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
//        return a;
//    }
//    vector<vec3> c = vector<vec3>(a.size());
//    for (unsigned int i=0; i < a.size(); i++)
//    {
//        c[i] = a[i] + b[i];
//    }
//    return c;
//}

//// Function for subtracting two std::vector<vec3>, with test for equal length
//vector<vec3> operator-(vector<vec3> a, vector<vec3> b){
//    if (a.size() != b.size())
//    {
//        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
//        return a;
//    }
//    vector<vec3> c = vector<vec3>(a.size());
//    for (unsigned int i=0; i < a.size(); i++)
//    {
//        c[i] = a[i] - b[i];
//    }
//    return c;
//}

*/

int main()
{
    // INPUT AREA

    // Set name of file to be read. File must contain 7 stats (mass, starting position (x,y,z),
    // starting velocity(vx,vy,vz)) on each line, separated by spaces. Each line counts as one body.
    const char* filename = "../ClusterData/solsyst.txt";

    // Galactic cluster(0) or Solar system(1)?
    int solsystOrCluster = 1;

    // Set integration method. RK4 is 1, adaptive RK4 is 2, Velocity Verlet is 3, Verlet is 4
    int method = 2;

    // Set endpoint of time calculations and timestep (dt)
    double total_time = 50;
    double timestep = 1e-1;

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
    cout << "Reading file " << filename << endl;
    cout << "Number of bodies = "<< astCluster.numberOfBodies()<<endl;

    // balle
    //        for(int i=0;i<astCluster.numberOfBodies();i++)
    //        {
    //            cout << astCluster.bodies.at(i).velocity << endl;
    //        }
    //const double pi = 4*std::atan(1.0); // Pi with double-precision
    //const double G = 4*pi*pi;

    //----------------------------------------------------------------------------------

    // FOR SOLAR SYSTEM ONLY:
    /* Correction for center of mass, allows input file
       to have the Sun at (0,0,0) with zero velocity and
       removes need to adjust center of mass when adding/removing
       bodies */

    if(solsystOrCluster==1)
    {
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

    else if(method == 4)  // Basic Störmer-Verlet
    {
        cout << "Chosen method is Basic Stormer-Verlet" << endl;
        Verlet::integrateVerlet(astCluster,timestep,n_steps);
    }

    // End timing
    finish = clock();
    double looptime = ((finish-start)/(double) CLOCKS_PER_SEC);

    cout <<endl<< "Simulation complete." << endl;
    cout <<endl<< "Time used  in loop: " << looptime << " s" << endl;
    return 0;
}


//balle
/* i matlab: dlmwrite

while global time<T
    do twice
        do twice
            update status of bodies using timestep dtmin
            global time ++
        update states of bodies using timestep dtmedium
    update states of bodies using timestep dtmax
    compute new timestep for all bodies


*/
