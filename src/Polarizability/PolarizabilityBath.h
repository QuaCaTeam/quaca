#ifndef POLARIZABILITYBATH_H
#define POLARIZABILITYBATH_H

#include <complex>
#include <cmath>
#include "MemoryKernel/MemoryKernelFactory.h"
#include "Polarizability.h"

class PolarizabilityBath : public Polarizability
{
private:
    // memory kernel needed to calculate alpha
    MemoryKernel *mu;

public:

    PolarizabilityBath(double omega_a, double alpha_zero, MemoryKernel *mu, GreensTensor *greens_tensor): Polarizability(omega_a, alpha_zero, greens_tensor) {this->mu = mu;};
    PolarizabilityBath(std::string input_file): Polarizability(input_file) {this->mu = MemoryKernelFactory::create(input_file);};

    std::complex<double> get_mu(double omega){return this->mu->mu(omega);};

    void calculate(cx_mat::fixed<3,3>& alpha, double omega);

};


#endif //POLARIZABILITYBATH_H
