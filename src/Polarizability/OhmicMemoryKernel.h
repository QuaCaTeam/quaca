#ifndef OHMICMEMORYKERNEL_H
#define OHMICMEMORYKERNEL_H

#include <complex>
#include <cmath>
#include "MemoryKernel.h"

class OhmicMemoryKernel : public MemoryKernel
{
  private:
    double gamma;

  public:
    std::complex<double> mu( double omega );
    OhmicMemoryKernel(double a);

};

#endif //OHMICMEMORYKERNEL_H
