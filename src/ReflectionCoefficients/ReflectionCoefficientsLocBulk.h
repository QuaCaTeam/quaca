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
  std::shared_ptr<Permittivity> permittivity;

public:
  /*!
   * Constructor for reflection coefficients of a local bulk medium.
   */
  ReflectionCoefficientsLocBulk(std::shared_ptr<Permittivity> permittivity);
  ReflectionCoefficientsLocBulk(const std::string &input_file);

  /*!
   * Returns the p- and s-polarized reflection coefficient.
   */
  void calculate(double omega, std::complex<double> kappa,
                 std::complex<double> &r_p,
                 std::complex<double> &r_s) const override;

  // getter functions
  std::complex<double> get_epsilon(double omega) const {
    return permittivity->calculate(omega);
  };

  std::complex<double> get_epsilon_omega(double omega) const {
    return permittivity->calculate_times_omega(omega);
  };

  // print info
  void print_info(std::ostream &stream) const override;
};

#endif // REFLECTIONCOEFFICIENTSLOCBULK_H
