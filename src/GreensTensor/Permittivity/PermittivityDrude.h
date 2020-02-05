#ifndef PERMITTIVITYDRUDE_H
#define PERMITTIVITYDRUDE_H

#include "Permittivity.h"
#include <complex>

//! A Drude model permittivity
class PermittivityDrude : public Permittivity {
private:
  double omega_p, gamma; // plasma frequency and damping coefficient

public:
  // constructors
  PermittivityDrude(double omega_p, double gamma);
  PermittivityDrude(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  /*!
   * Returns the numerical value of the permittivity scaled by omega.
   * @param Frequency
   */
  std::complex<double> epsilon_omega(double omega);

  /*!
   * Getter method for damping coefficient.
   */
  double get_gamma();

  /*!
   * Getter method for plasma frequency.
   */
  double get_omega_p();
};

#endif // PERMITTIVITYDRUDE_H
