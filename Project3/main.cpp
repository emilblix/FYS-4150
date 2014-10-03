#include <iostream>

using namespace std;

class Planet
{
protected:
    double Mass;        // Mass of planet in kg
    double Avg_dist;    // Average distance from the Sun in AU
public:
    double x;           // X-coordinate of planet in AU
    double y;           // Y-coordinate of planet in AU
    void setMass(double m)
    {
        Mass = m;
    }
    void setDist(double d)
    {
        Avg_dist = d;
    }
};

int main()
{
    Planet Earth;

    Earth.setMass(6e+24);
    Earth.setDist(1);

    int number_of_years = 2;        // Endpoint of time calculations
    int n_steps = 100;              // Number of calculation points
    double h = (double) number_of_years / (double) n_steps;

    cout<< h<<endl<<h*n_steps<<endl<<Earth.Mass<<endl;


    return 0;
}

