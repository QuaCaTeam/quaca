#include "MemoryKernel.h"


std::complex<double> OhmicMemoryKernel::mu(double omega)
{
    const   std::complex<double> I(0.0,1.0);  
    return gamma + I*0.0;
};
