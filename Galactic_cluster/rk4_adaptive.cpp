#include <rk4_adaptive.h>
#include <vec3.h>
#include <cmath>
#include <vector_operations.h>

using std::vector;

RK4_adaptive::RK4_adaptive()
{
}

void RK4_adaptive::RKF_45(Cluster &system, double dt_initial, double total_time)
{
    // Creating vectors
    int n_bodies = system.numberOfBodies();

    vector<vec3> A  = vector<vec3>(2*n_bodies);
    vector<vec3> K1 = vector<vec3>(2*n_bodies);
    vector<vec3> K2 = vector<vec3>(2*n_bodies);
    vector<vec3> K3 = vector<vec3>(2*n_bodies);
    vector<vec3> K4 = vector<vec3>(2*n_bodies);
    vector<vec3> K5 = vector<vec3>(2*n_bodies);
    vector<vec3> K6 = vector<vec3>(2*n_bodies);

    /*
    balle
    vector<int> test =vector<int>(5);
    for(int i=0;i<5;i++)
    {
        test[i]=i+6;
    }
    for (int i=0;i<test.size();i++)
    {
    std::cout << "test1    : "<<test[i]<<std::endl;
    }
    std::cout << "test1size: "<<test.size() << std::endl;
    test.clear();
    test.push_back(4);
    test.push_back(8);
    for (int i=0;i<test.size();i++)
    {
    std::cout << "test2    : "<<test[i] << std::endl;
    }
    std::cout << "test2size: "<<test.size() << std::endl;
    */

    for(int i=0;i<n_bodies;i++)
    {
        CelestialBody body = system.bodies[i];
        A[2*i]  = body.position;
        A[2*i+1]= body.velocity;
    }

//    std::cout << "A[2]: "<<A[2] << std::endl;
    double dt = dt_initial;
    double time_elapsed = 0;
    double tolerance = 1e-10;

    // Write first step to file for plotting
    system.dumpToFile(time_elapsed);

    // Time loop
    while(time_elapsed<=total_time)
    {
        // Print progress bar
        printf("Progress: %4.1f %% \r", 100.0*(time_elapsed)/(total_time));


        // Setting up A
        /* vector A = [position_1, velocity_1, position_2, velocity_2,  ... , position_(n_bodies), velocity_(n_bodies)]
                           (first body)            (second body)                 for all celestial bodies
        */
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody body = system.bodies[i];
            A[2*i]  = body.position;
            A[2*i+1]= body.velocity;
        }

        //----------------------------------------------------------------
        // Calculating K vectors, K1 is equal for RK4 and RKF45
        K1= dAdt(system,A);

        // REMOVE BEFORE DELIVERY? balle
        // Returning acceleration values to CelestialBody objects to track acceleration
        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            body->acceleration=K1[2*i+1];
        }

        // Calculating K vectors, continued
        // K1 = dt*dAdt(system,A)
        K1 = K1*dt;

        //---------------------------------------------------------------------------------------------
        // RK4 approximation using basic RK4 weighting

        K2 = dAdt(system, A+K1*(1/2.0)) * dt;

        K3 = dAdt(system, A+K2*(1/2.0)) * dt;

        K4 = dAdt(system, A+K3        ) * dt;

        vector<vec3> rk4approx = vector<vec3>(2*n_bodies);

        rk4approx = A + ((K1 + K2*2 + K3*2 + K4)*(1/6.0));

        //---------------------------------------------------------------------------------------------
        // RK5 approximation using RKF45 weighting

        K2 = dAdt(system, A + K1*(1/4.0) ) * dt;

        K3 = dAdt(system, A + K1*(3/32.0) + K2*(9/32.0)) * dt;

        K4 = dAdt(system, A + K1*(1932/2197.0) + K2*(-7200/2197.0) + K3*(7296/2197.0) ) * dt;

        K5 = dAdt(system, A + K1*(439/216.0) + K2*(-8) + K3*(3680/513.0) + K4*(-845/4104.0)) * dt;

        K6 = dAdt(system, A + K1*(-8/27.0) + K2*2 + K3*(-3544/2565.0) + K4*(1859/4104.0) + K5*(-11/40.0)) * dt;

        // RK5 approximation
        vector<vec3> rk5approx = vector<vec3>(2*n_bodies);

        rk5approx = A + K1*(16/135.0) + K3*(6656/12825.0) + K4*(28561/56430.0) + K5*(-9/50.0) + K6*(2/55.0);

        // RK4 approximation
        // Combining A + K1*25/216     + K3*1408/2565     + K4*2197/4101     - K5*1/5 into rk4approx
//         rk4approx =  A + K1*(25/216.0) + K3*(1408/2565.0) + K4*(2197/4101.0) + K5*(-1/5.0);



//        A = rk4approx;


        //-------------------------------------------------------------------------------------------------------
//        // Comparing total energy between methods
//        double currentEnergy = system.totalEnergy();
////        std::cout<<"currentEnergy: "<<currentEnergy<<std::endl;
////        currentEnergy = log10(fabs(currentEnergy));
////        std::cout<<"currentEnergy: "<<currentEnergy<<std::endl;

//        for(int i=0;i<n_bodies;i++)
//        {
//            CelestialBody *body = &system.bodies[i];
//            body->position = rk4approx[2*i];
//        }
////        double rk4energy = log10(fabs(system.totalEnergy()));
////        std::cout<<"rk4energy: "<<rk4energy<<std::endl;
//        double rk4energy =fabs(currentEnergy-rk4energy);
//        for(int i=0;i<n_bodies;i++)
//        {
//            CelestialBody *body = &system.bodies[i];
//            body->position = rk5approx[2*i];
//        }
////        double rk5energy = log10(fabs(system.totalEnergy()));
////        std::cout<<"rk5energy: "<<rk5energy<<std::endl;
//        double rk5energy = fabs(currentEnergy-rk5energy);

//        double diffMethods = fabs(rk4energy / rk5energy);
//        diffMethods = pow(10,diffMethods);
//        std::cout<<"diffm: "<<diffMethods<<std::endl;

        double diffMethods = 0;
        double rk5absvalue = 0;
        for (int i=0;i<2*n_bodies;i++)
        {
            vec3 diffVector = rk5approx[i]-rk4approx[i];
            diffMethods += fabs(diffVector.x()) + fabs(diffVector.y()) + fabs(diffVector.z());
            rk5absvalue += fabs(rk5approx[i].x()) + fabs(rk5approx[i].y()) + fabs(rk5approx[i].z());
        }
//        std::cout<<"rk5abs: "<<rk5absvalue<<std::endl;
        diffMethods = diffMethods/rk5absvalue;
//        std::cout<<"diffm: "<<diffMethods<<std::endl;
        diffMethods = diffMethods/(6*n_bodies);

        if(diffMethods > tolerance*10)
        {
            double s = pow(tolerance/(2*diffMethods),1/4.0);
//            s=0.5;
            dt = dt*s;
            printf("Time step decreased by a factor %4f to %4f \n",s,dt);
        }
        else if(diffMethods < tolerance*0.1)
        {
            double s = pow(tolerance/(2*diffMethods),1/4.0);
//            s=2;
            dt = dt*s;
            printf("Time step increased by a factor %4f to %4f \n",s,dt);
        }
        else
        {
            // Returning new position and velocity values to CelestialBody objects using the accepted RK5 approximation
            for(int i=0;i<n_bodies;i++)
            {
                CelestialBody *body = &system.bodies[i];
                body->position = rk5approx[2*i];
                body->velocity = rk5approx[2*i+1];
            }
            time_elapsed += dt;
            // Write to file for plotting
            system.dumpToFile(time_elapsed);
        }
        // The system needs to take a few steps at the start for the total energy to stabilize
//        else
//        {
////            A = A + rk5approx;
//            // Returning new position and velocity values to CelestialBody objects
//            for(int i=0;i<n_bodies;i++)
//            {
//                CelestialBody *body = &system.bodies[i];
//                body->position = rk5approx[2*i];
//                body->velocity = rk5approx[2*i+1];
//            }
//            time_elapsed += dt_initial;
//            // Write to file for plotting
//            system.dumpToFile(time_elapsed);
//        }
    }
}

vector<vec3> RK4_adaptive::dAdt(Cluster &system, vector<vec3> A)
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

/*
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

// Function for adding three std::vector<vec3>, with test for equal length
vector<vec3> RK4_adaptive::add3(vector<vec3> a, vector<vec3> b,vector <vec3> c)
{
    if (a.size() != b.size() || a.size() != c.size())
    {
        std::cout << "Error: Vectors a, b and c must be of equal length." << std::endl;
        return a;
    }
    vector<vec3> d = vector<vec3>(a.size());
    for (unsigned int i=0; i < a.size(); i++)
    {
        d[i] = a[i] + b[i] + c[i];
    }
    return d;
}

// Function for subtracting two std::vector<vec3>, with test for equal length
vector<vec3> RK4_adaptive::subtract(vector<vec3> a, vector<vec3> b)
{
    if (a.size() != b.size())
    {
        std::cout << "Error: Vectors a and b must be of equal length." << std::endl;
        return a;
    }
    vector<vec3> c = vector<vec3>(a.size());
    for (unsigned int i=0; i < a.size(); i++)
    {
        c[i] = a[i] - b[i];
    }
    return c;
}
*/


// balle
/*
vector<vec3> operator+( vector<vec3> &a, vector<vec3> &b)
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
*/
/*
#include <algorithm>
#include <functional>

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<T>());
    return result;
}
*/

/*
void RK4_adaptive::halfstep(Cluster &system, vector<vec3> A, vector<vec3> &K, double dt, int step)
{
    int testValue = 0; // Value to end loop when timestep is sufficiently low
    double newStepSize = dt;
    double minimumStepSize = 10e-6; // Lowest acceptable timestep value
    int halvingcounter = 0;         // Counts how many times the step size has been halved
    vector<vec3>   hsteppos;
    vector<double> hstepmass;
    vector<int>    hstepnum;

    // while loop to check
    while(testValue==0 && newStepSize > minimumStepSize)
    {
        hsteppos.clear();
        hstepmass.clear();
        hstepnum.clear();
        int n_bodies = system.numberOfBodies();

        for(int i=0;i<n_bodies;i++)
        {
            CelestialBody *body = &system.bodies[i];
            if(dt > 0.8/K[2*i+1].length())
            {
                hsteppos.push_back(K[2*i+1]);
                hstepmass.push_back(body->mass);
                hstepnum.push_back(i);
            }

        }
        if (hsteppos.size()>0.5)
        {
            newStepSize = newStepSize*0.5;
            step = step + newStepSize;
            halvingcounter += 1;


            // kalkulere dAdT
        }
        else {testValue=1;}
    }



}


//vector<vec3>RK4_adaptive::dHalfstepdt(Cluster system, vector<vec3> HS, vector<vec3> K1)
//{

//}

*/


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
