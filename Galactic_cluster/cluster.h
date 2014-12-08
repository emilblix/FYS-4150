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
    fstream outFile;

    Cluster();
    ~Cluster() {outFile.close();}
    void addCelestialBody(CelestialBody newBody);
    void calculateKineticAndPotentialEnergy();
    void calculateAcceleration();
    int numberOfBodies();
    double totalEnergy();
    void resetAllAcceleration();
    void dumpToFile(double time_elapsed);
};

#endif // CLUSTER_H
