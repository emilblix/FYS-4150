#include <rk_4.h>
#include <vec3.h>
#include <cmath>

using std::vector;

RK4::RK4()
{
}


//---------------------------------------------------------------------------------
// Using vector class vec3

void RK4::integrateCluster(Cluster &system, double dt, int n_steps)
{
    // Creating vectors
    int n_bodies = system.numberOfBodies();

    vector<vec3> A  = vector<vec3>(2*n_bodies);
    vector<vec3> K1 = vector<vec3>(2*n_bodies);
    vector<vec3> K2 = vector<vec3>(2*n_bodies);
    vector<vec3> K3 = vector<vec3>(2*n_bodies);
    vector<vec3> K4 = vector<vec3>(2*n_bodies);

    // Time loop
    for(int step=0;step<=n_steps;step++)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

        // Write to file for plotting
        system.dumpToFile(dt, step);

        // Setting up A
        /* vector A = [position_1, velocity_1,    (first body)
                   position_2, velocity_2,    (second body)
                   ... for all celestial bodies]
        */
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            A[2*i]  = body.position;
            A[2*i+1]= body.velocity;
        }

        // Runge-Kutta 4th order integration, start
        K1= dAdt(system,A);


        // REMOVE BEFORE DELIVERY? balle
        // Returning acceleration values to CelestialBody objects to track acceleration
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            body->acceleration=K1[2*i+1];
        }

        // Runge-Kutta 4th order integration, continued
        K1 = mult(K1,dt);
        K2 = dAdt(system,add(A,mult(K1,1/2.0))); K2 = mult(K2,dt);
        K3 = dAdt(system,add(A,mult(K2,1/2.0))); K3 = mult(K3,dt);
        K4 = dAdt(system,add(A,K3))            ; K4 = mult(K4,dt);

        // Combining (K1 + 2*K2 + 2*K3 + K4)/6 into K1
        K1 = add(K1,K4);
        K1 = add(K1,mult(K2,2));
        K1 = add(K1,mult(K3,2));
        K1 = mult(K1,1/6.0);

        A = add(A,K1);

        // Returning new position and velocity values to CelestialBody objects
        for(int i=1;i<n_bodies;i++)  // For stationary Sun, set startpoint i=1
        {
            CelestialBody *body = &system.bodies[i];
            body->position = A[2*i];
            body->velocity = A[2*i+1];
        }
    }
    // Print for last step
//        system.dumpToFile(dt,n_steps);
}


vector<vec3> RK4::dAdt(Cluster &system, vector<vec3> A)           // (24 * n(n+1)/2) + 3n= 12n^2 +15n FLOPs
{
    double pi = 4*std::atan(1.0);
    double G = 4*pi*pi;
    int n_bodies = system.numberOfBodies();
    vector<vec3> dAdt = vector<vec3>(2*n_bodies);

    // Zeroing the vector dAdt
    for(int i=0;i<2*n_bodies;i++)
    {
        dAdt[i].setToZero();
    }

    // Derivating vector A
    for(int i=0;i<n_bodies;i++)
    {
        // Derivative of position is velocity
        dAdt[2*i]   = A[2*i+1];

        double mass_i = system.bodies.at(i).mass;
        for(int j=i+1; j<n_bodies; j++)
        {
            // Calculating gravitational force between objects
            double mass_j = system.bodies.at(j).mass;
            vec3 dR = A[2*j] - A[2*i];          // dR pointing from body i to body j

            double dr = dR.length();
            double dr_cubed = dr*dr*dr;

            // For body i: a = mass_j/r^3*r
            double forcefactor = mass_j/dr_cubed;
            vec3 valueholder = dR*forcefactor;
            dAdt[2*i+1] = dAdt[2*i+1] + valueholder;

            // For body j: a = mass_i/r^3*r
            forcefactor = mass_i/dr_cubed;
            valueholder = dR*forcefactor;
            dAdt[2*j+1] = dAdt[2*j+1] - valueholder;
            //            std::cout<<"forcefactor = "<<i<<", "<<j<<": "<<valueholder<<std::endl;
        }
        // Multiplying with G = 4pi^2 at the end
        dAdt[2*i+1] = dAdt[2*i+1]*G;
    }
    return dAdt;
}

// Multiplication function for a std::vector<vec3> multiplied with a scalar
vector<vec3> RK4::mult(vector<vec3> a, double k)
{
    for (unsigned int i=0; i < a.size(); i++)
    {
        a[i] = a[i]*k;
    }
    return a;
}

// Function for adding two std::vector<vec3>, with test for equal length
vector<vec3> RK4::add(vector<vec3> a, vector<vec3> b)
{
    if (a.size() != b.size())
    {
        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
        return a;
    }
    vector<vec3> c = vector<vec3>(a.size());
    for (unsigned int i=0; i < a.size(); i++)
    {
        c[i] = a[i] + b[i];
    }
    return c;
}

/*
// Using x-, y- and z-coordinates as doubles

//void RK4::integrateCluster(Cluster &system, double dt)
//{
//    // Creating vectors
//    int n_bodies = system.numberOfBodies();

//    std::vector<double> A  = std::vector<double>(6*n_bodies);
//    std::vector<double> K1 = std::vector<double>(6*n_bodies);
//    std::vector<double> K2 = std::vector<double>(6*n_bodies);
//    std::vector<double> K3 = std::vector<double>(6*n_bodies);
//    std::vector<double> K4 = std::vector<double>(6*n_bodies);

//    // Setting up A
//     vector A = [x1, y1, z1, vx1, vy1, vz1    (first body)
//                   x2, y2, z2, vx2, vy2, vz2    (second body)
//                   ... for all celestial bodies
//
//    for(int i=0;i<n_bodies;i++)
//    {
//        CelestialBody body = system.bodies[i];
//        A[6*i]  = body.position.x();
//        A[6*i+1]= body.position.y();
//        A[6*i+2]= body.position.z();
//        A[6*i+3]= body.velocity.x();
//        A[6*i+4]= body.velocity.y();
//        A[6*i+5]= body.velocity.z();
//    }

//    // Runge-Kutta 4th order integration, start
//    K1= dAdt(system,A);

//    // Returning acceleration values to CelestialBody objects to track acceleration
//    for(int i=0;i<n_bodies;i++)
//    {
//        CelestialBody *body = &system.bodies[i];
//        body->acceleration.set(K1[6*i+3], K1[6*i+4], K1[6*i+5]);
//    }

//    // Runge-Kutta 4th order integration, continued
//    K1 = mult(K1,dt);
//    K2 = dAdt(system,add(A,mult(K1,1/2.0))); K2 = mult(K2,dt);
//    K3 = dAdt(system,add(A,mult(K2,1/2.0))); K3 = mult(K3,dt);
//    K4 = dAdt(system,add(A,K3))            ; K4 = mult(K4,dt);

//    // Combining (K1 + 2*K2 + 2*K3 + K4)/6 into K1
//    K1 = add(K1,K4);
//    K1 = add(K1,mult(K2,2));
//    K1 = add(K1,mult(K3,2));
//    K1 = mult(K1,1/6.0);

//    //    std::cout << "a before:" <<A[7]<<std::endl;

//    A = add(A,K1);

//    //std::cout<<K1[9]<<std::endl;
//    //    std::cout<< "a after:" << A[7]<<std::endl;


//    // Returning new position and velocity values to CelestialBody objects
//    for(int i=0;i<n_bodies;i++)  // For stationary Sun, set startpoint i=1
//    {
//        CelestialBody *body = &system.bodies[i];
//        body->position.set(A[6*i],   A[6*i+1], A[6*i+2]);
//        body->velocity.set(A[6*i+3], A[6*i+4], A[6*i+5]);
//    }
//}


//std::vector<double> RK4::dAdt(Cluster &system, std::vector<double> A)
//{
//    double pi = 4*std::atan(1.0);
//    double G = 4*pi*pi;
//    int n_bodies = system.numberOfBodies();
//    std::vector<double> dAdt = std::vector<double>(6*n_bodies);

//    // Zeroing the vector dAdt
//    for(int i=0;i<6*n_bodies;i++)
//    {
//        dAdt[i]=0;
//    }

//    // Derivating vector A
//    for(int i=0;i<n_bodies;i++)
//    {
//        // Derivative of position is velocity
//        dAdt[6*i]   = A[6*i+3];
//        dAdt[6*i+1] = A[6*i+4];
//        dAdt[6*i+2] = A[6*i+5];

//        double mass_i = system.bodies.at(i).mass;
//        for(int j=i+1; j<n_bodies; j++)
//        {
//            // Calculating gravitational force between objects
//            double mass_j = system.bodies.at(j).mass;
//            double dX = A[6*j]  -A[6*i];    // dX, dY and dZ "pointing" from body1 to body2
//            double dY = A[6*j+1]-A[6*i+1];
//            double dZ = A[6*j+2]-A[6*i+2];
//            double dr = std::sqrt(dX*dX+dY*dY+dZ*dZ);

//            double forcefactor = G*mass_i*mass_j/(dr*dr*dr);

//            dAdt[6*i+3] += dX*forcefactor/mass_i;
//            dAdt[6*i+4] += dY*forcefactor/mass_i;
//            dAdt[6*i+5] += dZ*forcefactor/mass_i;

//            dAdt[6*j+3] -= dX*forcefactor/mass_j;
//            dAdt[6*j+4] -= dY*forcefactor/mass_j;
//            dAdt[6*j+5] -= dZ*forcefactor/mass_j;
//        }
//    }
//    return dAdt;
//}

// Multiplication function for a std::vector<double> multiplied with a scalar
//std::vector<double> RK4::mult(std::vector<double> a, double k)
//{
//    for (unsigned int i=0; i < a.size(); i++) {
//        a[i] = a[i]*k;
//    }
//    return a;
//}


// Function for adding two std::vector<double>, with test for equal length
//std::vector<double> RK4::add(std::vector<double> a, std::vector<double> b)
//{

//    if (a.size() != b.size()) {
//        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
//        return a;
//    }
//    std::vector<double> c = std::vector<double>(a.size());
//    for (unsigned int i=0; i < a.size(); i++) {
//        c[i] = a[i] + b[i];
//    }
//    return c;
//}
*/
