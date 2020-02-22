#ifndef OHMICMEMORYKERNEL_H
#define OHMICMEMORYKERNEL_H

#include "MemoryKernel.h"
#include <cmath>
#include <complex>

//! An ohmic memory kernel
class OhmicMemoryKernel : public MemoryKernel {
private:
  double gamma; // damping coefficient

public:
  // constructors
  OhmicMemoryKernel(double gamma);
  OhmicMemoryKernel(std::string input_file, std::string section);
  OhmicMemoryKernel(std::string input_file);

  // calculate function
  std::complex<double> mu(double omega);

  // getter functions
  double get_gamma() { return this->gamma; };
};

#endif // OHMICMEMORYKERNEL_H
