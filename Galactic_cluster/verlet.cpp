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
    vector<vec3> lastAccel = vector<vec3>(n_bodies);
    vector<vec3> newAccel  = vector<vec3>(n_bodies);
    vector<vec3> position  = vector<vec3>(n_bodies);
    vector<vec3> velocity  = vector<vec3>(n_bodies);

    system.calculateAcceleration();
    for(int i=0;i<n_bodies;i++)
    {
        CelestialBody body = system.bodies[i];
        lastAccel[i] = body.acceleration;
        position[i]  = body.position;
        velocity[i]  = body.velocity;
    }

    // Write first step to file for plotting
    system.dumpToFile(0);

    // Time loop
    for(int step=1;step<=n_steps;step++)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

        // Calculate new position, split into two lines for readability
        position = add(position, mult(velocity,dt)         );     // r(i) = r(i-1) + v(i-1)*dt
        position = add(position, mult(lastAccel,dt*dt/2.0) );     // r(i) =    above line      + 1/2 a*dt^2

         // Returning new position values to CelestialBody objects
        for(int i=0;i<n_bodies;i++)  // For stationary Sun, set startpoint i=1
        {
            CelestialBody *body = &system.bodies[i];
            body->position = position[i];
        }
        // Calculate acceleration based on new positions
        system.calculateAcceleration();
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            newAccel[i] = body.acceleration;
        }
        // Calculate new velocity
        velocity = add(velocity, mult( add(lastAccel,newAccel), dt/2.0 )); // v(i) = v(i-1) + 1/2 [a(i-1)+a(i)]*dt

        // Returning new velocity values to CelestialBody objects
        for(int i=0;i<n_bodies;i++)  // For stationary Sun, set startpoint i=1
        {
            CelestialBody *body = &system.bodies[i];
            body->velocity = velocity[i];
        }

        // Pushing acceleration for use in next loop
        lastAccel = newAccel;

        // Write to file for plotting
        system.dumpToFile(dt*step);
    }

}

void Verlet::integrateVerlet(Cluster &system, double dt, int n_steps)
{

    // Setting initial vectors and variables
    int n_bodies = system.numberOfBodies();
    vector<vec3> lastPos      = vector<vec3>(n_bodies); // r(i-1)
    vector<vec3> currentPos   = vector<vec3>(n_bodies); // r(i)
    vector<vec3> newPos       = vector<vec3>(n_bodies); // r(i+1)
    vector<vec3> velocity     = vector<vec3>(n_bodies); // v(i)
    vector<vec3> acceleration = vector<vec3>(n_bodies); // a(i)

    system.calculateAcceleration();

    // Write first step to file for plotting
    system.dumpToFile(0);

    for(int i=0;i<n_bodies;i++)
    {
        CelestialBody body = system.bodies[i];
        lastPos[i]      = body.position;
        velocity[i]     = body.velocity;
        acceleration[i] = body.acceleration;
    }
    // Calculate position at step 1, split into two lines for readability
    currentPos = add(lastPos,mult(velocity,dt));                // r(1) = r(0) + v(0)dt
    currentPos = add(currentPos,mult(acceleration,dt*dt/2.0));  // r(1) =  above line   + 1/2 a(0) dt^2


    // Time loop
    for(int step=1;step<=n_steps;step++)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

        // Calculate acceleration based on new positions
        system.calculateAcceleration();
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            acceleration[i] = body.acceleration;
        }

        // Calculate new position, split into two lines for readability
        newPos = subtract(mult(currentPos,2), lastPos); // r(i+1) = 2r(i) - r(i-1)
        newPos = add(newPos, mult(acceleration,dt*dt)); // r(i+1) =  above line    + a(i)*dt^2


        // Calculate new velocity
        velocity = mult(subtract(newPos,lastPos), 1/(2*dt)); // v(i) = (r(i+1) - r(i-1))/2dt

        //Returning new position and velocity values to CelestialBody objects
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            body->position = newPos[i];
            body->velocity = velocity[i];
        }

        // Write to file for plotting
        system.dumpToFile(dt*step);

        // Pushing positions for use in next loop
        lastPos = currentPos;
        currentPos = newPos;

    }

}

// Multiplication function for a std::vector<vec3> multiplied with a scalar
vector<vec3> Verlet::mult(vector<vec3> a, double k)
{
    for (unsigned int i=0; i < a.size(); i++)
    {
        a[i] = a[i]*k;
    }
    return a;
}

// Function for adding two std::vector<vec3>, with test for equal length
vector<vec3> Verlet::add(vector<vec3> a, vector<vec3> b)
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

// Function for subtracting two std::vector<vec3>, with test for equal length
vector<vec3> Verlet::subtract(vector<vec3> a, vector<vec3> b)
{
    if (a.size() != b.size())
    {
        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
        return a;
    }
    vector<vec3> c = vector<vec3>(a.size());
    for (unsigned int i=0; i < a.size(); i++)
    {
        c[i] = a[i] - b[i];
    }
    return c;
}
