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
  OhmicMemoryKernel(const std::string& input_file, const std::string& section);
  OhmicMemoryKernel(const std::string& input_file);

  // calculate function
  std::complex<double> mu(double omega);

  // getter functions
  double get_gamma() const { return this->gamma; };

  void print_info(std::ofstream &file);
};

#endif // OHMICMEMORYKERNEL_H
