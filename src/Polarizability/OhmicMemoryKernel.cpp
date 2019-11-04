#include "OhmicMemoryKernel.h"

OhmicMemoryKernel::OhmicMemoryKernel(double a)
{
  this->gamma = a;
};

std::complex<double> OhmicMemoryKernel::mu(double omega)
{
    const   std::complex<double> gammac(this->gamma,0E0);  
    return gammac;
};
