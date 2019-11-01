#ifndef POLARIZABILITY_H
#define POLARIZABILITY_H

#include <complex>
#include <cmath>

#include "MemoryKernel.h"

class Polarizability
{
private:

    // parameters
    double omega_a;
    double alpha_zero;

    // memory kernel and greens tensor needed to calculate alpha
    MemoryKernel *mu;

    // array where alpha is actually stored
    std::complex<double> alpha[3][3];

public:

    Polarizability(double a, double b, MemoryKernel *mu);


};


#endif //POLARIZABILITY_H
