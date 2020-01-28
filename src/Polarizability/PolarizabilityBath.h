#ifndef POLARIZABILITYBATH_H
#define POLARIZABILITYBATH_H

#include "MemoryKernel/MemoryKernelFactory.h"
#include "Polarizability.h"
#include <cmath>
#include <complex>

class PolarizabilityBath : public Polarizability {
private:
  MemoryKernel *mu; // memory kernel needed to calculate alpha

public:
  // constructors
  PolarizabilityBath(double omega_a, double alpha_zero, MemoryKernel *mu,
                     GreensTensor *greens_tensor);
  PolarizabilityBath(std::string input_file);

  // calculate the polarizability tensor
  void calculate_tensor(cx_mat::fixed<3, 3> &alpha,
                        Options_Polarizability opts);

  // getter function for memory kernel
  std::complex<double> get_mu(double omega) { return this->mu->mu(omega); };
};

#endif // POLARIZABILITYBATH_H
