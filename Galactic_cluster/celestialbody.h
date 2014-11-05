#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H
#include <vec3.h>

class CelestialBody
{

public:
    double mass;                        // Relative mass of body compared to the Sun
    vec3 position;                      // Coordinates of body in AU
    vec3 velocity;                      // Velocity of body in AU/yr
    vec3 acceleration;                  // Acceleration of body in AU/yr²
    vec3 force;                         // Force on body in dimensionless variable
    CelestialBody(double mass_, double x, double y, double z, double vx, double vy, double vz);
    CelestialBody(vec3 pos, vec3 vel, double mass_);
    void resetForce();
    void resetAcceleration();
};

#endif // CELESTIAL_BODY_H
