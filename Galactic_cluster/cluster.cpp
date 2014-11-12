#include <cluster.h>
#include <iomanip>
#include <celestialbody.h>
#include <cmath>

Cluster::Cluster()
{
    outFile.open("./Data/pos.dat", ios::out);
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();
}

void Cluster::addCelestialBody(CelestialBody newBody)
{
    bodies.push_back(newBody);
}

void Cluster::calculateKineticAndPotentialEnergy() // Calculates kinetic and potential energy of system
{
    const double pi = 4*std::atan(1.0); // Pi with double-precision
    const double G = 4*pi*pi;
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();


    for(int i=0; i<numberOfBodies(); i++)
    {
        CelestialBody body1 = bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++)
        {
            CelestialBody body2 = bodies[j];

            // Distance between bodies
            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();

            // Add potential energy to total
            potentialEnergy -= G*body1.mass*body2.mass/dr;
        }

        kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
    }
}

void Cluster::calculateForces() // Calculates forces between bodies
{
    // Set forces to zero to avoid addition from previous
    resetAllForces();
    for(int i=0; i<numberOfBodies(); i++)
    {
        CelestialBody *body1 = &bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++)
        {
            CelestialBody *body2 = &bodies[j];
            vec3 deltaRVector = body2->position - body1->position;    // deltaRVector pointing from body1 to body2
            double dr = deltaRVector.length();

            double forcefactor = body1->mass*body2->mass/(dr*dr*dr);
            vec3 forceVector= deltaRVector*forcefactor;

            body1->force = body1->force+forceVector;     // Force on body1 points same direction as deltaRVector
            body2->force = body2->force-forceVector;     // Force on body2 points opposite direction as deltaRvector
        }
    }
}


void Cluster::calculateAcceleration() // Calculates forces between bodies
{
    const double pi = 4*std::atan(1.0); // Pi with double-precision
    const double G = 4*pi*pi;
    // Set forces to zero to avoid addition from previous run
    resetAllAcceleration();
    for(int i=0; i<numberOfBodies(); i++)
    {
        CelestialBody *body1 = &bodies.at(i);
        for(int j=i+1; j<numberOfBodies(); j++)
        {
            CelestialBody *body2 = &bodies[j];
            vec3 deltaRVector = body2->position - body1->position;    // deltaRVector pointing from body1 to body2
            double dr = deltaRVector.length();
            double dr_cubed = dr*dr*dr;

            double accFactor = body2->mass/dr_cubed;
            vec3 accVector= deltaRVector*accFactor;
            body1->acceleration = body1->acceleration+accVector;     // Force on body1 points same direction as deltaRVector

            accFactor = body1->mass/dr_cubed;
            accVector= deltaRVector*accFactor;
            body2->acceleration = body2->acceleration-accVector;     // Force on body2 points opposite direction as deltaRvector
        }

        // Finally, multiplying the sum of "acceleration vectors" with 4pi^2 to get the proper acceleration
        body1->acceleration = body1->acceleration*G;
    }
}

int Cluster::numberOfBodies() // Gives the number of bodies in the system
{
    return bodies.size();
}

double Cluster::totalEnergy()       // Returns sum of kinetic and potential energy in system
{
    calculateKineticAndPotentialEnergy();
    return kineticEnergy + potentialEnergy;
}

void Cluster::resetAllForces()      //  Sets the force vector of all celestial bodies to zero
{
    for(int i=0; i<numberOfBodies(); i++)
    {
        bodies[i].resetForce();
    }
}

void Cluster::resetAllAcceleration()      //  Sets the force vector of all celestial bodies to zero
{
    for(int i=0; i<numberOfBodies(); i++)
    {
        bodies[i].resetAcceleration();
    }
}

vec3 Cluster::forceAtPosition(int bodyNumber, vec3 pos)
{
    CelestialBody body1 = bodies[bodyNumber];
    vec3 forceOnBody = vec3(0,0,0) ;
    for(int i=0; i<numberOfBodies(); i++)
    {
        if(i != bodyNumber)
        {
            CelestialBody body2 = bodies[i];
            vec3 deltaRVector = body2.position - pos;          // deltaRVector pointing from forceBody to body2
            double dr = deltaRVector.length();                  // Distance between bodies
            double forcefactor = body1.mass*body2.mass/(dr*dr*dr);
            vec3 forceVector= deltaRVector*forcefactor;
            forceOnBody = forceOnBody+forceVector;              // Force on forceBody points same direction as deltaRVector
        }

        // Calculate the potential energy here
    }
    return forceOnBody;
}

void Cluster::dumpToFile(double timestep, int step) // Updates the file "pos.dat" with x and y positions of all bodies
{
    outFile << timestep*step << " " << totalEnergy() << " ";
    for (int i = 0 ;i<numberOfBodies(); i++)
    {
        CelestialBody body = bodies[i];
        outFile << std::setprecision(16)<<body.position.x() << " " << body.position.y() << " ";

//        outFile << body.position.x() << " " << body.position.y() << " " << body.velocity.x() << " " << body.velocity.y() << " "<< body.acceleration.x() << " " << body.acceleration.y() << " ";
    }
    outFile << "\n";
}

