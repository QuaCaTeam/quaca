#ifndef POLARIZABILITY_H
#define POLARIZABILITY_H

#include <complex>
#include <cmath>
#include <armadillo>

#include "MemoryKernel/MemoryKernel.h"

using namespace arma;

class Polarizability
{
protected:
    // parameters
    double omega_a;
    double alpha_zero;

    // array where alpha is actually stored
    cx_mat::fixed<3,3> alpha;

public:

    virtual cx_mat::fixed<3,3> calculate(double omega) =0;

    // getter methods
    double get_omega_a();
    double get_alpha_zero();
};


#endif //POLARIZABILITY_H
