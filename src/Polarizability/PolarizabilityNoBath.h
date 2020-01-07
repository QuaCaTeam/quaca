#ifndef POLARIZABILITYNOBATH_H
#define POLARIZABILITYNOBATH_H

#include <complex>
#include <cmath>
#include "Polarizability.h"

class PolarizabilityNoBath : public Polarizability
{
public:

    PolarizabilityNoBath(double omega_a, double alpha_zero, GreensTensor *greens_tensor): Polarizability(omega_a, alpha_zero, greens_tensor) {};
    PolarizabilityNoBath(std::string input_file): Polarizability(input_file) {};

    void calculate(cx_mat::fixed<3,3>& alpha, double omega);

};


#endif //POLARIZABILITYNOBATH_H
