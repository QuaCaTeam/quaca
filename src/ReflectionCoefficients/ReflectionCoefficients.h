#ifndef REFLECTIONCOEFFICIENTS_H
#define REFLECTIONCOEFFICIENTS_H

#include <armadillo>
#include <cmath>
#include <complex>

// abstract class for reflection coefficients
class ReflectionCoefficients {
public:
  // returns the reflection coefficients
  virtual void ref(std::complex<double> &r_p, std::complex<double> &r_s,
                   double omega, std::complex<double> kappa) = 0;

  virtual void print_info(std::ofstream &file) = 0;
};

#endif // REFLECTIONCOEFFICIENTS_H
