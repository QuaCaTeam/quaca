#ifndef REFLECTIONCOEFFICIENTS_H
#define REFLECTIONCOEFFICIENTS_H

#include <armadillo>
#include <cmath>
#include <complex>
//! An abstract permittivity class
/*!
 * This is an abstract class for reflection coefficients
 */
class ReflectionCoefficients {
public:
  /*!
   * Returns the reflection coefficients
   * @param omega Frequency
   */
  virtual void ref(std::complex<double> &r_p, std::complex<double> &r_s, double omega, std::complex<double> kappa) = 0;
};

#endif // REFLECTIONCOEFFICIENTS_H
