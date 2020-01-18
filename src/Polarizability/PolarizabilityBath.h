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

    PolarizabilityBath(double omega_a, double alpha_zero, MemoryKernel *mu, GreensTensor *greens_tensor);
    PolarizabilityBath(std::string input_file);

    void calculate_tensor(cx_mat::fixed<3,3>& alpha, Options_Polarizability opts);

    std::complex<double> get_mu(double omega){return this->mu->mu(omega);};

};


#endif //POLARIZABILITYBATH_H
