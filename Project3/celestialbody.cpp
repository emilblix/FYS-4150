#include <celestialbody.h>
#include <vec3.h>



CelestialBody::CelestialBody(double x, double y, double vx, double vy, double mass_, double radius_,char const *name_)
{
    position = vec3(x,y,0);
    velocity = vec3(vx,vy,0);
    mass = mass_;
    name = name_;
    radius = radius_;
}

CelestialBody::CelestialBody(vec3 pos, vec3 vel, double mass_)
{
    position = pos;
    velocity = vel;
    mass = mass_;
    //resetForce();
}
void CelestialBody::resetForce()
{
    force.setToZero();
}
