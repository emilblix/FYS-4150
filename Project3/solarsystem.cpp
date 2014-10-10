#include <solarsystem.h>
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
//            std::cout<<1<<forceVector.x()<<std::endl;
            body1->force = body1->force-forceVector;     // Force on body1 points opposite direction as deltaRVector
            body2->force = body2->force+forceVector;     // Force on body2 points same direction as deltaRvector
//            std::cout<<2<<body2->force.x()<<std::endl;
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
void SolarSystem::dumpToFile(double timestep, int step) // Updates the file "pos.dat" with x and y positions of all bodies
{
    outFile << timestep*step << " ";
    for (int i = 0 ; i < numberOfBodies(); i++)
    {
        CelestialBody *body = bodies[i];
        outFile << body->position.x() << " " << body->position.y() << " ";
    }
    outFile << "\n";
}

