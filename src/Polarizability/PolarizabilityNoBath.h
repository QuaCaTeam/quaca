#ifndef POLARIZABILITYNOBATH_H
#define POLARIZABILITYNOBATH_H

#include "Polarizability.h"
#include <cmath>
#include <complex>
#include <memory>

class PolarizabilityNoBath : public Polarizability {
public:
  PolarizabilityNoBath(double omega_a, double alpha_zero,
                       std::shared_ptr<GreensTensor> greens_tensor);
  PolarizabilityNoBath(const std::string &input_file);

  void calculate_tensor(double omega, cx_mat::fixed<3, 3> &alpha,
                        Tensor_Options fancy_complex) const;
};

#endif // POLARIZABILITYNOBATH_H
