#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include <celestialbody.h>
#include <vector>
#include <fstream>
#include <iostream>

using std::vector;
using std::fstream;
using std::ios;

class SolarSystem
{
public:
    vector<CelestialBody*> bodies;
    double kineticEnergy;
    double potentialEnergy;
    vec3 angularMomentum;
    fstream outFile;

    SolarSystem();
    ~SolarSystem() {outFile.close();}
    void addCelestialBody(CelestialBody &newBody);
    void calculateForcesAndEnergy();
    int numberOfBodies();
    double totalEnergy();
    void resetAllForces();
    void dumpToFile(double timestep, int step);
};

#endif // SOLARSYSTEM_H
