#include <runge_kutta_4.h>

RK4::RK4()
{
}

void RK4::integrateSolarSystem(SolarSystem &system, double dt)
{
    // Do RK4 integration here
    int n_bodies = system.numberOfBodies();
    system.calculateForcesAndEnergy();

    for(int i=1;i<n_bodies;i++)
    {
//        char *variableName = new char[1000];          // Lage egen class for RK-vektorer?
//        sprintf(variableName, "K1_%01d", i);
//        vec3 variableName;
//        sprintf(variableName, "K2_%01d", i);
//        vec3 variableName;
//        sprintf(variableName, "K3_%01d", i);
//        vec3 variableName;
//        sprintf(variableName, "K4_%01d", i);
//        vec3 variableName;
        CelestialBody *body = system.bodies[i];
//        vec3 K1,K2,K3,K4;

        K1_1=system.forceAtPosition(body,body->position);
//        Lag funksjon force(position) eller acceleration(position) og funksjon
    std::cout<<"K1 = "<<K1_1<<endl;
    }



//    vec RK4(vec A, double dt){
//        vec K1,K2,K3,K4;
//        K1 = dt * force(A);
//        K2 = dt * force(A+K1/2);
//        K3 = dt * force(A+K2/2);
//        K4 = dt * force(A+K3);

//        return A +(K1 + 2*K2 + 2*K3 + K4)/6;
}
