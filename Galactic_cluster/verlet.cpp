#include <verlet.h>
#include <vec3.h>
#include <cmath>
#include <vector_operations.h>

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

        // Calculate new position
        // r(i)  =  r(i-1)  +   v(i-1)*dt +         a* dt^2*1/2
        position = position + velocity*dt + lastAccel*(dt*dt/2.0);

         // Returning new position values to CelestialBody objects
        for(int i=0;i<n_bodies;i++)
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
        // v(i)  =  v(i-1)  + [  a(i-1) +  a(i)  ]* dt*1/2
        velocity = velocity + (lastAccel+newAccel)*(dt/2.0);

        // Returning new velocity values to CelestialBody objects
        for(int i=0;i<n_bodies;i++)
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
    // Calculate position at step 1
    // r(1)    =   r(0)  +     v(0)*dt +     a(0)    * dt^2*1/2
    currentPos = lastPos + velocity*dt + acceleration*(dt*dt/2.0);

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

        // Calculate new position
        // r(i+1) =    r(i)*2   - r(i-1)  +    a(i)     * dt^2
        newPos = currentPos*2.0 - lastPos + acceleration*(dt*dt);

        // Calculate new velocity
        // v(i) =  (r(i+1)- r(i-1))  /   (2*dt)
        velocity = (newPos-lastPos) * (1/(2*dt));

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
