#include <celestialbody.h>
#include <vec3.h>



CelestialBody::CelestialBody(double mass_, double x, double y, double z, double vx, double vy, double vz)
{
    position = vec3(x,y,z);
    velocity = vec3(vx,vy,vz);
    mass = mass_;
}

CelestialBody::CelestialBody(vec3 pos, vec3 vel, double mass_)
{
    position = pos;
    velocity = vel;
    mass = mass_;
}

void CelestialBody::resetAcceleration()
{
    acceleration.setToZero();
}
