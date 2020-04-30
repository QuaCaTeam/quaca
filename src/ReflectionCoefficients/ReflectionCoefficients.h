#ifndef REFLECTIONCOEFFICIENTS_H
#define REFLECTIONCOEFFICIENTS_H

#include <armadillo>
#include <cmath>
#include <complex>

// abstract class for reflection coefficients
class ReflectionCoefficients {
public:
  // returns the reflection coefficients
  virtual void calculate(double omega, std::complex<double> kappa,
                         std::complex<double> &r_p,
                         std::complex<double> &r_s) const = 0;
};

#endif // REFLECTIONCOEFFICIENTS_H
