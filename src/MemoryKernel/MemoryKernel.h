#ifndef MEMORYKERNEL_H
#define MEMORYKERNEL_H

#include <cmath>
#include <complex>

//! An abstract class for memory kernels
class MemoryKernel {
public:
  // Returns the memory kernel given a frequency omega.
  virtual std::complex<double> mu(double omega) = 0;

  virtual void print_info(std::ofstream &file) = 0;
};

#endif // MEMORYKERNEL_H
