#ifndef POLARIZABILITYBATH_H
#define POLARIZABILITYBATH_H

#include <complex>
#include <cmath>
#include "Polarizability.h"

class PolarizabilityBath : public Polarizability
{
private:
    // memory kernel and greens tensor needed to calculate alpha
    MemoryKernel *memorykernel;

public:

    PolarizabilityBath(double a, double b, MemoryKernel *mu);
    PolarizabilityBath(std::string input_file);

    std::complex<double> get_mu(double omega);

    void calculate(cx_mat::fixed<3,3>& alpha, double omega);

};


#endif //POLARIZABILITYBATH_H
