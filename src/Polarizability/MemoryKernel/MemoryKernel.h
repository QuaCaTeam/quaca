#ifndef MEMORYKERNEL_H
#define MEMORYKERNEL_H

#include <cmath>
#include <complex>

//! An abstract class for memory kernels
class MemoryKernel {
protected:
  std::string type; // type of memory kernel

public:
  // constructor
  MemoryKernel();
  MemoryKernel(std::string input_file);

  // Returns the memory kernel given a frequency omega.
  virtual std::complex<double> mu(double omega) = 0;
};

#endif // MEMORYKERNEL_H
