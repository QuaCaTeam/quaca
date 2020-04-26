#ifndef MEMORYKERNEL_H
#define MEMORYKERNEL_H

#include <cmath>
#include <complex>

//! An abstract class for memory kernels
class MemoryKernel {
public:
  // Returns the memory kernel given a frequency omega.
  virtual std::complex<double> calculate(double omega) = 0;
};

#endif // MEMORYKERNEL_H
