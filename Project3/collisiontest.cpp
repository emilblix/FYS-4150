#include <collisiontest.h>
#include <cmath>


int collisionTest(SolarSystem &system)
{
    // Measures the distance between all bodies, prints all colliding bodies and breaks main loop in case of collision
    int collisiontest=0;
    for (int i=0;i<system.numberOfBodies();i++)
    {
        CelestialBody *body1 = system.bodies[i];
        for (int j=i+1;j<system.numberOfBodies();j++)
        {
            CelestialBody *body2 = system.bodies[j];
            vec3 deltaRVector = body1->position - body2->position;    // deltaRVector pointing from body2 to body1
            double dr = deltaRVector.length();
            if(fabs(dr)<2*(body1->radius + body2->radius))          // Collision if distance between objects is less than twice their combined radius
            {
                if(collisiontest==1)
                {
                    printf(" between %s and %s,",body1->name, body2->name);
                }
                else
                {
                    collisiontest=1;
                    printf("Collision between %s and %s,",body1->name, body2->name);
                }
            }
        }
    }
    return collisiontest;
}
