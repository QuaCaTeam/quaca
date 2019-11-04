#ifndef MEMORYKERNEL_H
#define MEMORYKERNEL_H

#include <complex>
#include <cmath>

class MemoryKernel
{
public:
    virtual std::complex<double> mu( double omega ) =0;
};

#endif //MEMORYKERNEL_H
