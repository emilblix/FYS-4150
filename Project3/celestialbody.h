#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H
#include <vec3.h>

class CelestialBody
{

public:
    char const* name;
    double mass;                        // Relative mass of body compared to the Sun
    double radius;                      // Radius in AU
    vec3 position;                      // Coordinates of planet in AU
    vec3 velocity;                      // Velocity of planet in AU/yr
    vec3 acceleration;                  // Acceleration of body in AU/yr²
    vec3 force;                         // Force on planet in dimensionless variable
    CelestialBody(double x, double y, double vx, double vy, double mass_, double radius_, const char *name_);
    CelestialBody(vec3 pos, vec3 vel, double mass_);
    void resetForce();
};

#endif // CELESTIAL_BODY_H
