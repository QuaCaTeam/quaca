#ifndef POLARIZABILITYBATH_H
#define POLARIZABILITYBATH_H

#include "../MemoryKernel/MemoryKernelFactory.h"
#include "Polarizability.h"
#include <cmath>
#include <complex>

class PolarizabilityBath : public Polarizability {
private:
  // memory kernel needed to calculate alpha
  std::shared_ptr<MemoryKernel> mu;

public:
  PolarizabilityBath(double omega_a, double alpha_zero,
                     std::shared_ptr<MemoryKernel> mu,
                     std::shared_ptr<GreensTensor> greens_tensor);
  PolarizabilityBath(const std::string &input_file);

  void calculate_tensor(double omega, cx_mat::fixed<3, 3> &alpha,
                        Tensor_Options fancy_complex) const;

  // getter function for memory kernel
  std::complex<double> get_mu(double omega) const {
    return mu->calculate(omega);
  };
};

#endif // POLARIZABILITYBATH_H
