#ifndef REFLECTIONCOEFFICIENTSLOCSLAB_H
#define REFLECTIONCOEFFICIENTSLOCSLAB_H

#include "../Permittivity/PermittivityFactory.h"
#include "ReflectionCoefficients.h"
#include <armadillo>
#include <complex>

//! Reflection coefficients for a local bulk medium
/*!
 * This is a class implementing a Drude model permittivity, which is given by
 * \f$ \varepsilon(\omega) = 1 - \frac{\omega_p^2}{\omega (\omega + i
 * \gamma)}\f$ Its attributes are the damping coefficient and the plasma
 * frequency
 */
class ReflectionCoefficientsLocSlab : public ReflectionCoefficients {
private:
  // permittivity is needed to describe the surface's response
  Permittivity *permittivity;
  double thickness;

public:
  /*!
   * Constructor for reflection coefficients of a local bulk medium.
   */
  ReflectionCoefficientsLocSlab(Permittivity *permittivity, double thickness);

  ReflectionCoefficientsLocSlab(std::string input_file);

  /*!
   * Returns the p- and s-polarized reflection coefficient.
   */
  void ref(std::complex<double> &r_p, std::complex<double> &r_s, double omega, std::complex<double> kappa);
  // getter functions
  std::complex<double> get_epsilon(double omega) {
    return this->permittivity->epsilon(omega);
  };
  double get_thickness(){
   return this->thickness;
  };
  };

#endif // REFLECTIONCOEFFICIENTSLOCSLAB_H
