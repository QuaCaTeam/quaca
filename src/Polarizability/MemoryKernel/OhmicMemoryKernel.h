#ifndef OHMICMEMORYKERNEL_H
#define OHMICMEMORYKERNEL_H

#include <complex>
#include <cmath>
#include "MemoryKernel.h"

//! An ohmic memory kernel
/*!
* This is a class implementing an ohmic memory kernel, i.e. a memory kernel
* given by \f$ \mu(\omega) = \gamma \f$.
* It takes the damping coefficient as an input
*/
class OhmicMemoryKernel : public MemoryKernel
{
private:
  double gamma;

public:

  // constructors
  OhmicMemoryKernel(double gamma): gamma(gamma) {};
  OhmicMemoryKernel(std::string input_file);

  // getter functions
  std::complex<double> mu( double omega );
  double get_gamma(){return this->gamma;};

};

#endif //OHMICMEMORYKERNEL_H
