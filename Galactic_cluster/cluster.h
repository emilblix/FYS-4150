#ifndef CLUSTER_H
#define CLUSTER_H

#include <celestialbody.h>
#include <vector>
#include <fstream>
#include <iostream>

using std::vector;
using std::fstream;
using std::ios;

class Cluster
{
public:
    vector<CelestialBody> bodies;
    double kineticEnergy;
    double potentialEnergy;
    vec3 angularMomentum;
    fstream outFile;

    Cluster();
    ~Cluster() {outFile.close();}
    void addCelestialBody(CelestialBody newBody);
    void calculateKineticAndPotentialEnergy();
    void calculateForces();
    int numberOfBodies();
    double totalEnergy();
    void resetAllForces();
    vec3 forceAtPosition(int bodyNumber, vec3 pos);
    void dumpToFile(double timestep, int step);
};

#endif // CLUSTER_H
