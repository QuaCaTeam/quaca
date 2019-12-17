#ifndef INTEGRATIONS_H
#define INTEGRATIONS_H
#include <iostream>
#include <cmath>
#include <armadillo>

/*!
* Wrapper function for the cquad integration routine of the gsl.
* @param my_f Function to integrate (must have form of gsl_function)
* @param a lower integration bound
* @param b upper integration bound
* @param relerr relative error
* @param relerr absolute error
*/
double cquad(double my_f(double, void *), void* params, double a , double b, double relerr, double epsabs);

/*!
* Wrapper function for the qags integration routine of the gsl.
* @param my_f Function to integrate (must have form of gsl_function)
* @param a lower integration bound
* @param b upper integration bound
* @param relerr relative error
* @param relerr absolute error
*/
double qags(double my_f(double, void *), void* params, double a , double b, double relerr, double epsabs);

/*!
* Wrapper function for the qagiu integration routine of the gsl.
* @param my_f Function to integrate (must have form of gsl_function)
* @param a lower integration bound
* @param relerr relative error
* @param relerr absolute error
*/
double qagiu(double my_f(double, void *), void* params, double a , double relerr, double epsabs);

#endif //INTEGRATIONS_H
