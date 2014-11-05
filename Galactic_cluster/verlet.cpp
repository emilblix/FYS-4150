#include <verlet.h>
#include <vec3.h>
#include <cmath>

using std::vector;

Verlet::Verlet()
{
}


void Verlet::velocityVerlet(Cluster &system, double dt, int n_steps)
{

    // Setting initial vectors and variables
    int n_bodies = system.numberOfBodies();
    std::vector<vec3> lastAccel = std::vector<vec3>(n_bodies);
    std::vector<vec3> newAccel  = std::vector<vec3>(n_bodies);
    std::vector<vec3> position  = std::vector<vec3>(n_bodies);
    std::vector<vec3> velocity  = std::vector<vec3>(n_bodies);

    system.calculateAcceleration();
    for(int i=0;i<n_bodies;i++)
    {
        CelestialBody body = system.bodies[i];
        lastAccel[i] = body.acceleration;
        position[i]  = body.position;
        velocity[i]  = body.velocity;
    }

    // Write first step to file for plotting
    system.dumpToFile(dt, 0);

    // Time loop
    for(int step=1;step<=n_steps;step++)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

        // Calculate new position, split into two lines for readability
        position = add(position, mult(velocity,dt)         );     // r(i) = r(i-1) + v(i-1)*dt
        position = add(position, mult(lastAccel,dt*dt/2.0) );     // r(i) =    above line      + 1/2 a*dt^2

        // Calculate acceleration based on new positions
        system.calculateAcceleration();
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            newAccel[i] = body.acceleration;
        }
        // Calculate new velocity
        velocity = add(velocity, mult( add(lastAccel,newAccel), dt/2.0 )); // v(i) = v(i-1) + 1/2 [a(i-1)+a(i)]*dt

        // Returning new position and velocity values to CelestialBody objects
        for(int i=0;i<n_bodies;i++)  // For stationary Sun, set startpoint i=1
        {
            CelestialBody *body = &system.bodies[i];
            body->position = position[i];
            body->velocity = velocity[i];
        }

        // Pushing acceleration for use in next loop
        lastAccel = newAccel;

        // Write to file for plotting
        system.dumpToFile(dt, step);
    }

}

void Verlet::integrateVerlet(Cluster &system, double dt, int n_steps)
{

    // Setting initial vectors and variables
    int n_bodies = system.numberOfBodies();
    std::vector<vec3> lastPos      = std::vector<vec3>(n_bodies);
    std::vector<vec3> newPos       = std::vector<vec3>(n_bodies);
    std::vector<vec3> currentPos   = std::vector<vec3>(n_bodies);
    std::vector<vec3> velocity     = std::vector<vec3>(n_bodies);
    std::vector<vec3> acceleration = std::vector<vec3>(n_bodies);

    system.calculateAcceleration();
    for(int i=0;i<n_bodies;i++)
    {
        CelestialBody body = system.bodies[i];
        lastPos[i] = body.position;
        currentPos[i].setToZero();
        velocity[i].setToZero();
    }

    // Write first step to file for plotting
    system.dumpToFile(dt, 0);

    // Time loop
    for(int step=1;step<=n_steps;step++)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

        // Calculate new position, split into two lines for readability
        position = add(position, mult(velocity,dt)         );     // r(i) = r(i-1) + v(i-1)*dt
        position = add(position, mult(lastAccel,dt*dt/2.0) );     // r(i) =    above line      + 1/2 a*dt^2

        // Calculate acceleration based on new positions
        system.calculateAcceleration();
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            newAccel[i] = body.acceleration;
        }
        // Calculate new velocity
        velocity = add(velocity, mult( add(lastAccel,newAccel), dt/2.0 )); // v(i) = v(i-1) + 1/2 [a(i-1)+a(i)]*dt

        // Returning new position and velocity values to CelestialBody objects
        for(int i=0;i<n_bodies;i++)  // For stationary Sun, set startpoint i=1
        {
            CelestialBody *body = &system.bodies[i];
            body->position = position[i];
            body->velocity = velocity[i];
        }

        // Pushing acceleration for use in next loop
        lastAccel = newAccel;

        // Write to file for plotting
        system.dumpToFile(dt, step);
    }

}

// Multiplication function for a std::vector<vec3> multiplied with a scalar
std::vector<vec3> Verlet::mult(std::vector<vec3> a, double k)
{
    for (unsigned int i=0; i < a.size(); i++) {
        a[i] = a[i]*k;
    }
    return a;
}

// Function for adding two std::vector<vec3>, with test for equal length
std::vector<vec3> Verlet::add(std::vector<vec3> a, std::vector<vec3> b)
{

    if (a.size() != b.size()) {
        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
        return a;
    }
    std::vector<vec3> c = std::vector<vec3>(a.size());
    for (unsigned int i=0; i < a.size(); i++) {
        c[i] = a[i] + b[i];
    }
    return c;
}
