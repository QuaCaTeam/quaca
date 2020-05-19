#ifndef MEMORYKERNEL_H
#define MEMORYKERNEL_H

#include <cmath>
#include <complex>

//! An abstract class for memory kernels
class MemoryKernel {
public:
  // Returns the memory kernel given a frequency omega.
  virtual std::complex<double> calculate(double omega) const = 0;

  // print info
  virtual void print_info(std::ostream &stream) const =0;
};

#endif // MEMORYKERNEL_H
