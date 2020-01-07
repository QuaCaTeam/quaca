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
  PermittivityDrude(double omega_p, double gamma): omega_p(omega_p), gamma(gamma) {};
  PermittivityDrude(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // getter functions
  double get_gamma() const{return this->gamma;};
  double get_omega_p()const {return this->omega_p;};
};

#endif // PERMITTIVITYDRUDE_H
