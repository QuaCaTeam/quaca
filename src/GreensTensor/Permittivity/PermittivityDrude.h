#ifndef PERMITTIVITYDRUDE_H
#define PERMITTIVITYDRUDE_H

#include <complex>
#include "Permittivity.h"

//! A Drude model permittivity
/*!
* This is a class implementing a Drude model permittivity, which is given by
* \f$ \varepsilon(\omega) = 1 - \frac{\omega_p^2}{\omega (\omega + i \gamma)}\f$
* Its attributes are the damping coefficient and the plasma frequency
*/
class PermittivityDrude : public Permittivity
{
private:
  double omega_p, gamma;

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
