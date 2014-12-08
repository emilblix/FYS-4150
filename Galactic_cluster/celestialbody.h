#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H
#include <vec3.h>

class CelestialBody
{

public:
    double mass;                        // Relative mass of body compared to the Sun
    vec3 position;                      // Coordinates of body
    vec3 velocity;                      // Velocity of body
    vec3 acceleration;                  // Acceleration of body
    CelestialBody(double mass_, double x, double y, double z, double vx, double vy, double vz);
    CelestialBody(vec3 pos, vec3 vel, double mass_);
    void resetAcceleration();
};

#endif // CELESTIAL_BODY_H
