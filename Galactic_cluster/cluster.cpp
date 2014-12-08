#include <cluster.h>
#include <iomanip>
#include <celestialbody.h>
#include <cmath>


Cluster::Cluster()
{
    outFile.open("../ClusterData/pos.dat", ios::out);
    kineticEnergy = 0;
    potentialEnergy = 0;
}

void Cluster::addCelestialBody(CelestialBody newBody)
{
    bodies.push_back(newBody);
}

void Cluster::calculateKineticAndPotentialEnergy() // Calculates kinetic and potential energy of system
{
    double G = gravitationalConstant;
    kineticEnergy = 0;
    potentialEnergy = 0;

    for(int i=0; i<numberOfBodies(); i++)
    {
        CelestialBody body1 = bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++)
        {
            CelestialBody body2 = bodies[j];

            // Distance between bodies
            vec3 deltaRVector = body1.position - body2.position;
            double dr = 1e-6+deltaRVector.length();

            // Add potential energy to total
            potentialEnergy -= G*body1.mass*body2.mass/dr;
        }

        kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
    }
}

void Cluster::calculateAcceleration() // Calculates forces between bodies
{
    double G = gravitationalConstant;
    // Set acceleration to zero to avoid addition from previous run
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
        // Finally, multiplying the sum of "acceleration vectors" with G to get the proper acceleration
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

void Cluster::resetAllAcceleration()      //  Sets the force vector of all celestial bodies to zero
{
    for(int i=0; i<numberOfBodies(); i++)
    {
        bodies[i].resetAcceleration();
    }
}

void Cluster::dumpToFile(double time_elapsed) // Updates the file "pos.dat" with x and y positions of all bodies
{
    calculateKineticAndPotentialEnergy();
    outFile << time_elapsed << " "<< kineticEnergy << " " << potentialEnergy << " ";
    for (int i = 0 ;i<numberOfBodies(); i++)
    {
        CelestialBody body = bodies[i];
        outFile << std::setprecision(16)<<body.position.x() << " " << body.position.y() << " "<< body.position.z()<<" ";
    }
    outFile << "\n";
}

