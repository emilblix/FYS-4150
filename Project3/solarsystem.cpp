#include <solarsystem.h>
#include <iomanip>
#include <celestialbody.h>
#include <cmath>

SolarSystem::SolarSystem()
{
    outFile.open("pos.dat", ios::out);
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();
}

void SolarSystem::addCelestialBody(CelestialBody &newBody)
{
    bodies.push_back(&newBody);
}

void SolarSystem::calculateKineticAndPotentialEnergy() // Calculates kinetic and potential energy of system
{
    const double pi = 4*std::atan(1.0); // Pi with double-precision
    const double G = 4*pi*pi;
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();


    for(int i=0; i<numberOfBodies(); i++)
    {
        CelestialBody *body1 = bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++)
        {
            CelestialBody *body2 = bodies[j];

            // Distance between bodies
            vec3 deltaRVector = body1->position - body2->position;
            double dr = deltaRVector.length();

            // Add potential energy to total
            potentialEnergy -= G*body1->mass*body2->mass/dr;
        }

        kineticEnergy += 0.5*body1->mass*body1->velocity.lengthSquared();
    }
}

void SolarSystem::calculateForces() // Calculates forces between bodies
{
    // Set forces to zero to avoid addition from previous
    resetAllForces();
    for(int i=0; i<numberOfBodies(); i++)
    {
        CelestialBody *body1 = bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++)
        {
            CelestialBody *body2 = bodies[j];
            vec3 deltaRVector = body2->position - body1->position;    // deltaRVector pointing from body1 to body2
            double dr = deltaRVector.length();

            double forcefactor = body1->mass*body2->mass/(dr*dr*dr);
            vec3 forceVector= deltaRVector*forcefactor;

            body1->force = body1->force+forceVector;     // Force on body1 points same direction as deltaRVector
            body2->force = body2->force-forceVector;     // Force on body2 points opposite direction as deltaRvector
        }
    }
}

int SolarSystem::numberOfBodies() // Gives the number of bodies in the solar system
{
    return bodies.size();
}

double SolarSystem::totalEnergy()       // Returns sum of kinetic and potential energy in solar system
{
    calculateKineticAndPotentialEnergy();
    return kineticEnergy + potentialEnergy;
}
void SolarSystem::resetAllForces()      //  Sets the force vector of all celestial bodies to zero
{
    for(int i=0; i<numberOfBodies(); i++)
    {
        bodies[i]->resetForce();
    }
}

vec3 SolarSystem::forceAtPosition(int bodyNumber, vec3 pos)
{
    CelestialBody *body1 = bodies[bodyNumber];
    vec3 forceOnBody = vec3(0,0,0) ;
    for(int i=0; i<numberOfBodies(); i++)
    {
        if(i != bodyNumber)
        {
            CelestialBody *body2 = bodies[i];
            vec3 deltaRVector = body2->position - pos;          // deltaRVector pointing from forceBody to body2
            double dr = deltaRVector.length();                  // Distance between bodies
            double forcefactor = body1->mass*body2->mass/(dr*dr*dr);
            vec3 forceVector= deltaRVector*forcefactor;
            forceOnBody = forceOnBody+forceVector;              // Force on forceBody points same direction as deltaRVector
        }

        // Calculate the potential energy here
    }
    return forceOnBody;
}

void SolarSystem::dumpToFile(double timestep, int step) // Updates the file "pos.dat" with x and y positions of all bodies
{
    outFile << timestep*step << " " << totalEnergy() << " ";
    for (int i = 0 ; i < numberOfBodies(); i++)
    {
        CelestialBody *body = bodies[i];
        outFile << std::setprecision(16)<<body->position.x() << " " << body->position.y() << " ";

//        outFile << body->position.x() << " " << body->position.y() << " " << body->velocity.x() << " " << body->velocity.y() << " "<< body->acceleration.x() << " " << body->acceleration.y() << " ";
    }
    outFile << "\n";
}

