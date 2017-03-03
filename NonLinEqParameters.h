//file NonLinEqParameters.h
#include <iostream>
#include <string>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <errno.h>

using namespace std;

class non_lin_eq_parameters                              //members by default are private
{
      typedef double (*pfn)(double,double *); //define a function pointer to a one variable function with parameters contains in an array
      double delta;                         //root tolerance, return when updated interval is narrower than it
      double epsn;                         //residual tolerance, return when residual is smaller than it
      int    maxit;                        //maximum number of iterations allowed
      double *varargin;

public:                                     //public members
       non_lin_eq_parameters(double a,double b, int c,double *param, int size)   //constructor
       {

            varargin = new double[size];
            delta = a; epsn = b; maxit = c;
            for (int r = 0; r < size; r++ )
            {
                varargin[r]=param[r];
            }

       }
       ~non_lin_eq_parameters() {delete[] varargin;}              //destructor

       double secant(double ,double  ,pfn);
       double BrentDekker(double,double,pfn);
};


double non_lin_eq_parameters::secant(double x0,double x1 ,pfn f)
{
/*******************************************************************************
secant's algm:  find an approximate root of f(x) = 0
x0 :            an initial guess of a root
x1 :            second initial guess of a root
f:              the function whose root is to be found
delta:          program stops when distance of two iterates is < it
epsn:           program stops when residual is less than it
maxit:          maximum number of iterations allowed
*******************************************************************************/

      double v0 = f(x0,varargin);
      double v1 = f(x1,varargin);                           //fcn value at initial guess

      double xnew;                               //store new iterate

      for (int k = 1; k <= maxit; k++)
      {
          double derv = (v1-v0)/(x1-x0);                 //derivative at xp
          if (!derv)
          {
                    cout << "division by 0 occurred in secant()." << endl;
                    exit(1); //stop if divisor == 0
          }

          xnew = x1 - v1/derv; //compute new iterate
          v0 = v1;
          v1 = f(xnew,varargin); //fcn value at iterate
          if (fabs(xnew-x1)<delta || fabs(v1)<epsn) return xnew;
          x0 = x1;
          x1 = xnew;
      }
      return xnew;
}

double non_lin_eq_parameters::BrentDekker(double x0,double x1,pfn f)
{
/*******************************************************************************
Brent-Dekker's algm:  find an approximate root of f(x) = 0
x0 :            an initial guess of a root
x1 :            second initial guess of a root
f:              the function whose root is to be found
delta:          program stops when distance of two iterates is < it
epsn:           program stops when residual is less than it
maxit:          maximum number of iterations allowed
*******************************************************************************/

 double a;
 double b;
 double c;
 double d;
 double e;
 double m;
 double p;
 double q;
 double r;
 double s;

 a  = x0;
 b  = x1;

 double fa;
 double fb;
 double fc;

 fa = f(a,varargin);
 fb = f(b,varargin);


 double tol;

// if there is no sign change between f(a) and f(b), we are searchning for one
 if (fa*fb > 0)
 {

    double fp1;
    double fp2;
    for (int u = 1; u <= 1000; u++)
    {
        fp1 = f(a+0.1*u,varargin);

        if (fp1*fb < 0)
        {
           fa = fp1;
           a +=0.1*u;
           break;
        }
        fp2 = f(a-0.1*u,varargin);

        if (fp2*fb < 0)
        {
           fa = fp2;
           a-= 0.1*u;
           break;
        }
        if (u == 1000)
        {
              cout << "The function must changes sign in the initial interval" << endl;
              cout << "The function does not change with the initial point a varied up to 1000%" << endl;
              exit(1);
        }
    }


 }

 c  = a;
 fc = fa;
 d  = b - c;
 e  = d;



 for (int k = 1; k <= maxit; k++)
 {

     // The three points a, b et c must satisfy:
     //    f(a)*f(b) < 0
     //    abs(f(b)) <= abs(f(a))
     //    c = previous b (we could have c = a)
     //    The nest point will be either:
     //    a bissection approx at (a+b)/2
     //    a secant approx calculated with b and c
     //    an approximation of the inverse quadratic interpolation calculated with a, b and c

     if (fa*fb > 0)
     {
        a = c;     fa = fc;
        d = b - c; e  = d;
     }

     if ( fabs(fa) < fabs(fb))
     {
        c = b;   b  = a;  a  = c;
        fc = fb; fb = fa; fa = fc;
     }


    //Stopping criterion

    m = 0.5*(a - b);
    tol = 2.0*delta*max(fabs(b),1.0);

    if ((fabs(m) <= tol) || (fabs(fb) <= epsn))return b;

    //We choose between bissection, secant or IQI
    if ((fabs(e) < tol) || (fabs(fc) <= fabs(fb)))
    {
       // Bissection
       d = m;
       e = m;
    }
    else
    {
        s = fb/fc;
        if (a == c)
        {
           // Secant
           p = 2.0*m*s;
           q = 1.0 - s;
        }
        else
        {
            //IQI
            q = fc/fa;
            r = fb/fa;
            p = s*(2.0*m*q*(q - r) - (b - c)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        }

        if (p > 0)
            q = -q;
        else
            p = -p;

        //Est-ce que le point de secante ou d'IQI est acceptable?
        //Do we accept the point we got from the secant or IQI
        if ((2.0*p < 3.0*m*q - fabs(tol*q)) && (p < fabs(0.5*e*q)))
        {
           e = d;
           d = p/q;
        }
        else
        {
            d = m;
            e = m;
        }
    }

    //Next point
    c  = b;
    fc = fb;

    if (fabs(d) > tol)
       b = b + d;
    else
        b = b - (b-a)/fabs(b-a)*tol;

    fb = f(b,varargin);


  }
  return b;
}
