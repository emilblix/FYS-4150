#include <collisiontest.h>
#include <cmath>


int collisionTest(SolarSystem &system)
{
    int collisiontest=0;
    for (unsigned int i=0;i<system.numberOfBodies();i++)
    {
        CelestialBody *body1 = system.bodies[i];
        for (unsigned int j=i+1;j<system.numberOfBodies();j++)
        {
            CelestialBody *body2 = system.bodies[j];
            vec3 deltaRVector = body1->position - body2->position;    // deltaRVector pointing from body2 to body1
            double dr = deltaRVector.length();
            if(fabs(dr)<0.005)
            {
                collisiontest=1;
                printf("Collision between %s and %s,",body1->name, body2->name);
                break;
            }
        }
        if(collisiontest==1)
        {
            break;
        }
    }
    return collisiontest;
}
