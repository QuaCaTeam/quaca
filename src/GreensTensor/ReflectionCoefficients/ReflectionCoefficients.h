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
  virtual void refvec(cx_vec ref(2), double omega) = 0;
};

#endif // REFLECTIONCOEFFICIENTS_H
