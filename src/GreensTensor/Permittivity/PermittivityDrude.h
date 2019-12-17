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
  /*! plasma frequency */
  double omega_p;

  /*! damping coefficient */
  double gamma;

public:

  /*!
  * Constructor for a Drude model permittivity.
  * @param a damping coefficient
  * @param b plasma frequency
  */
  PermittivityDrude(double a, double b);

  /*!
  * Constructor for a Drude model permittivity fron an input file.
  * @param input_file .ini file containing the parameters
  */
  PermittivityDrude(std::string input_file);

  /*!
  * Returns the numerical value of the permittivity.
  * @param Frequency
  */
  std::complex<double> epsilon(double omega);

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
