#ifndef REFLECTIONCOEFFICIENTSLOCBULK_H
#define REFLECTIONCOEFFICIENTSLOCBULK_H

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
class ReflectionCoefficientsLocBulk : public ReflectionCoefficients {
private:
  // permittivity is needed to describe the surface's response
  Permittivity *permittivity;

public:
  /*!
   * Constructor for reflection coefficients of a local bulk medium.
   */
  ReflectionCoefficientsLocBulk(Permittivity *permittivity);
  ReflectionCoefficientsLocBulk(const std::string &input_file);

  /*!
   * Returns the p- and s-polarized reflection coefficient.
   */
  void ref(std::complex<double> &r_p, std::complex<double> &r_s, double omega,
           std::complex<double> kappa);

  // getter functions
  std::complex<double> get_epsilon(double omega) const {
    return this->permittivity->epsilon(omega);
  };


  void print_info(std::ofstream &file);
};

#endif // REFLECTIONCOEFFICIENTSLOCBULK_H
