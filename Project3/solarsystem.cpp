#include <solarsystem.h>
#include <iomanip>
#include <celestialbody.h>

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

void SolarSystem::calculateForcesAndEnergy() // Calculates forces between bodies
{
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();


    for(int i=0; i<numberOfBodies(); i++)
    {
        CelestialBody *body1 = bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++)
        {
            CelestialBody *body2 = bodies[j];
            vec3 deltaRVector = body1->position - body2->position;    // deltaRVector pointing from body2 to body1
            double dr = deltaRVector.length();

            double forcefactor = body1->mass*body2->mass/(dr*dr*dr);
            vec3 forceVector= deltaRVector*forcefactor;

            body1->force = body1->force-forceVector;     // Force on body1 points opposite direction as deltaRVector
            body2->force = body2->force+forceVector;     // Force on body2 points same direction as deltaRvector

            // Calculate the potential energy here
        }

        kineticEnergy += 0.5*body1->mass*body1->velocity.lengthSquared();
    }
}

int SolarSystem::numberOfBodies() // Gives the number of bodies in the solar system
{
    return bodies.size();
}

double SolarSystem::totalEnergy()       // Returns sum of kinetic and potential energy in solar system
{
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
    outFile << timestep*step << " ";
    for (int i = 1 ; i < numberOfBodies(); i++) // Startpoint 0 to include Sun
    {
        CelestialBody *body = bodies[i];
//        outFile << std::setprecision(16)<<body->position.x() << " " << body->position.y() << " ";

        outFile << body->position.x() << " " << body->position.y() << " " << body->velocity.x() << " " << body->velocity.y()
        << " "<< body->acceleration.x() << " " << body->acceleration.y() << " ";
    }
    outFile << "\n";
}

