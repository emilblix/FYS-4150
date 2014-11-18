#include <rk4_adaptive.h>
#include <vec3.h>
#include <cmath>

using std::vector;

RK4_adaptive::RK4_adaptive()
{
}

void RK4_adaptive::integrateClusterAdaptive(Cluster &system, double dt, int n_steps)
{
    // Creating vectors
    int n_bodies = system.numberOfBodies();

    vector<vec3> A  = vector<vec3>(2*n_bodies);
    vector<vec3> K1 = vector<vec3>(2*n_bodies);
    vector<vec3> K2 = vector<vec3>(2*n_bodies);
    vector<vec3> K3 = vector<vec3>(2*n_bodies);
    vector<vec3> K4 = vector<vec3>(2*n_bodies);


//    for(int i=0;i<n_bodies;i++)
//    {
//        CelestialBody body = system.bodies[i];
//        A[2*i]  = body.position;
//        A[2*i+1]= body.velocity;
//    }

    std::cout << "A[6]: "<<A[6] << std::endl;

    // Time loop
    for(int step=0;step<=n_steps;step++)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*((double)step)/((double)n_steps));

        // Write to file for plotting
        system.dumpToFile(dt, step);

        // Setting up A
        /* vector A = [position_1, velocity_1,    (first body)
                   position_2, velocity_2,    (second body)
                   ... for all celestial bodies]
        */
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            A[2*i]  = body.position;
            A[2*i+1]= body.velocity;
        }

        // Runge-Kutta 4th order integration, start
        K1= dAdt(system,A);

        // REMOVE BEFORE DELIVERY? balle
        // Returning acceleration values to CelestialBody objects to track acceleration
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            body->acceleration=K1[2*i+1];
        }


        // Runge-Kutta 4th order integration, continued
        K1 = mult(K1,dt);
        K2 = dAdt(system,add(A,mult(K1,1/2.0))); // dAdt(system, A+1/2*dt*K1)
        K2 = mult(K2,dt);
        K3 = dAdt(system,add(A,mult(K2,1/2.0))); // dAdt(system, A+1/2*dt*K2)
        K3 = mult(K3,dt);
        K4 = dAdt(system,add(A,K3))            ; // dAdt(system, A+dt*K3)

        RK4_adaptive::calcHalfstep(system,K4,dt,step);

        K4 = mult(K4,dt);

        // Combining (K1 + 2*K2 + 2*K3 + K4)/6 into K1
        K1 = add(K1,K4);
        K1 = add(K1,mult(K2,2));
        K1 = add(K1,mult(K3,2));
        K1 = mult(K1,1/6.0);

        A = add(A,K1);

        // Returning new position and velocity values to CelestialBody objects
        for(int i=1;i<n_bodies;i++)  // For stationary Sun, set startpoint i=1
        {
            CelestialBody *body = &system.bodies[i];
            body->position = A[2*i];
            body->velocity = A[2*i+1];
        }
    }
}

void RK4_adaptive::calcHalfstep(Cluster &system, vector<vec3> &K4, double dt, int step)
{
    int testValue = 0; // Value to end loop when timestep is sufficiently low
    double newStepSize = dt;
    double minimumStepSize = 10e-6; // Lowest acceptable timestep value
    do{
        vector<vec3>   halfstep;
        vector<double> hstepmass;
        vector<int>    hsteppos;
        int n_bodies = system.numberOfBodies();

        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            if(dt > 0.8/K4[2*i+1].length())
            {
                halfstep.push_back(K4[2*i+1]);
                hstepmass.push_back(body->mass);
                hsteppos.push_back(i);
            }

        }
        if (halfstep.size()>1.0)
        {
            newStepSize = newStepSize*0.5;
            step = step + newStepSize;


            // kalkulere dAdT
        }
        else {testValue=1;}
    }while(testValue==0 && newStepSize > minimumStepSize);



}


//vector<vec3>RK4_adaptive::dHalfstepdt(Cluster system, vector<vec3> HS, vector<vec3> K1)
//{

//}



vector<vec3> RK4_adaptive::dAdt(Cluster &system, vector<vec3> A)           // (24 * n(n+1)/2) + 3n= 12n^2 +15n FLOPs
{
    double pi = 4*std::atan(1.0);
    double G = 4*pi*pi;
    int n_bodies = system.numberOfBodies();
    vector<vec3> dAdt = vector<vec3>(2*n_bodies);

    // Zeroing the vector dAdt
    for(int i=0;i<2*n_bodies;i++)
    {
        dAdt[i].setToZero();
    }

    // Derivating vector A
    for(int i=0;i<n_bodies;i++)
    {
        // Derivative of position is velocity
        dAdt[2*i]   = A[2*i+1];

        double mass_i = system.bodies.at(i).mass;
        for(int j=i+1; j<n_bodies; j++)
        {
            // Calculating gravitational force between objects
            double mass_j = system.bodies.at(j).mass;
            vec3 dR = A[2*j] - A[2*i];          // dR pointing from body i to body j

            double dr = dR.length();
            double dr_cubed = dr*dr*dr;

            // For body i: a = mass_j/r^3*r
            double forcefactor = mass_j/dr_cubed;
            vec3 valueholder = dR*forcefactor;
            dAdt[2*i+1] = dAdt[2*i+1] + valueholder;

            // For body j: a = mass_i/r^3*r
            forcefactor = mass_i/dr_cubed;
            valueholder = dR*forcefactor;
            dAdt[2*j+1] = dAdt[2*j+1] - valueholder;
        }
        // Multiplying with G = 4pi^2 at the end
        dAdt[2*i+1] = dAdt[2*i+1]*G;
    }
    return dAdt;
}

// Multiplication function for a std::vector<vec3> multiplied with a scalar
vector<vec3> RK4_adaptive::mult(vector<vec3> a, double k)
{
    for (unsigned int i=0; i < a.size(); i++)
    {
        a[i] = a[i]*k;
    }
    return a;
}

// Function for adding two std::vector<vec3>, with test for equal length
vector<vec3> RK4_adaptive::add(vector<vec3> a, vector<vec3> b)
{
    if (a.size() != b.size())
    {
        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
        return a;
    }
    vector<vec3> c = vector<vec3>(a.size());
    for (unsigned int i=0; i < a.size(); i++)
    {
        c[i] = a[i] + b[i];
    }
    return c;
}


// rk4_adaptive.cpp

/*
  Function to advance set of coupled first-order o.d.e.s by single step
  using adaptive fourth-order Runge-Kutta scheme

     x       ... independent variable
     y       ... array of dependent variables
     h       ... step-length
     t_err   ... actual truncation error per step
     acc     ... desired truncation error per step
     S       ... step-length cannot change by more than this factor from
                  step to step
     rept    ... number of step recalculations
     maxrept ... maximum allowable number of step recalculations
     h_min   ... minimum allowable step-length
     h_max   ... maximum allowable step-length
     flag    ... controls manner in which truncation error is calculated

  Requires right-hand side routine

        void rhs_eval (double x, Array<double,1> y, Array<double,1>& dydx)

  which evaluates derivatives of y (w.r.t. x) in array dydx.

  Function advances equations by single step whilst attempting to maintain
  constant truncation error per step of acc:

    flag = 0 ... error is absolute
    flag = 1 ... error is relative
    flag = 2 ... error is mixed

  If step-length falls below h_min then routine aborts
*/
/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <blitz/array.h>

using namespace blitz;

void rk4_fixed (double&, Array<double,1>&,
                void (*)(double, Array<double,1>, Array<double,1>&),
                double);

void rk4_adaptive (double& x, Array<double,1>& y,
                   void (*rhs_eval)(double, Array<double,1>, Array<double,1>&),
                   double& h, double& t_err, double acc,
                   double S, int& rept, int maxrept,
                   double h_min, double h_max, int flag)
{
  // Array y assumed to be of extent n,  where n is no. of coupled
  // equations
  int n = y.extent(0);

  // Declare local arrays
  Array<double,1> y0(n), y1(n);

  // Declare repetition counter
  static int count = 0;

  // Save initial data
  double x0 = x;
  y0 = y;

  // Take full step
  rk4_fixed (x, y, rhs_eval, h);

  // Save data
  y1 = y;

  // Restore initial data
  x = x0;
  y = y0;

  // Take two half-steps
  rk4_fixed (x, y, rhs_eval, h/2.);
  rk4_fixed (x, y, rhs_eval, h/2.);

  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error
      for (int i = 0; i < n; i++)
        {
          err = fabs (y(i) - y1(i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < n; i++)
        {
          err = fabs ((y(i) - y1(i)) / y(i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else
    {
      // Use mixed truncation error
      for (int i = 0; i < n; i++)
        {
          err1 = fabs ((y(i) - y1(i)) / y(i));
          err2 = fabs (y(i) - y1(i));
          err = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err == 0.) {t_err = 1.e-15;}

  // Calculate new step-length
  double h_est = h * pow (fabs (acc / t_err), 0.2);

  // Prevent step-length from changing by more than factor S
  if (h_est / h > S)
    h *= S;
  else if (h_est / h < 1. / S)
    h /= S;
  else
    h = h_est;

  // Prevent step-length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;

  // Abort if step-length falls below h_min
  if (fabs(h) < h_min)
    {
      printf ("Error - |h| < hmin\n");
      exit (1);
    }

  // If truncation error acceptable take step
  if ((t_err <= acc) || (count >= maxrept))
    {
      rept = count;
      count = 0;
    }
  // If truncation error unacceptable repeat step
  else
    {
      count++;
      x = x0;
      y = y0;
      rk4_adaptive (x, y, rhs_eval, h, t_err, acc,
                    S, rept, maxrept, h_min, h_max, flag);
    }

  return;
}
*/
