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
    double gravitationalConstant;
    vec3 angularMomentum;
    fstream outFile;

    Cluster();
    ~Cluster() {outFile.close();}
    void addCelestialBody(CelestialBody newBody);
    void calculateKineticAndPotentialEnergy();
    void calculateForces();
    void calculateAcceleration();
    int numberOfBodies();
    double totalEnergy();
    void resetAllForces();
    void resetAllAcceleration();
    vec3 forceAtPosition(int bodyNumber, vec3 pos);
    void dumpToFile(double time_elapsed);
};

#endif // CLUSTER_H
