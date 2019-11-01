#ifndef MEMORYKERNEL_H
#define MEMORYKERNEL_H

#include <complex>
#include <cmath>

class MemoryKernel
{
public:
    virtual std::complex<double> mu( double omega ) =0;
};

class OhmicMemoryKernel : public MemoryKernel
{
private:
    double gamma;

public:
    std::complex<double> mu( double omega );

};

#endif //MEMORYKERNEL_H
