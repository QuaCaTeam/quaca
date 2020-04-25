#ifndef INTEGRATIONS_H
#define INTEGRATIONS_H
#include <armadillo>
#include <cmath>
#include <iostream>

// wrapper functions for integration routines of the gsl
double cquad(double my_f(double, void *), void *params, double a, double b,
             double relerr, double epsabs);
double qags(double my_f(double, void *), void *params, double a, double b,
            double relerr, double epsabs);
double qagiu(double my_f(double, void *), void *params, double a, double relerr,
             double epsabs);

#endif // INTEGRATIONS_H
