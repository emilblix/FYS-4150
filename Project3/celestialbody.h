#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H
#include <vec3.h>

class CelestialBody
{

public:
    char* name;
    double mass;                        // Relative mass of body compared to the Sun
    double Avg_dist;                    // Average distance from the Sun in AU
    vec3 position;                      // Coordinates of planet in AU
    vec3 velocity;                      // Velocity of planet in AU/yr
    vec3 acceleration;                  // Acceleration of body in AU/yr²
    vec3 force;                         // Force on planet in dimensionless variable
    CelestialBody(vec3 pos, vec3 vel, double mass_);
    CelestialBody(double x, double y, double vx, double vy, double mass_, char* name_);
    void resetForce();
};

#endif // CELESTIAL_BODY_H
