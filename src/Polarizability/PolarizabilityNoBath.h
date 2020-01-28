#ifndef POLARIZABILITYNOBATH_H
#define POLARIZABILITYNOBATH_H

#include "Polarizability.h"
#include <cmath>
#include <complex>

class PolarizabilityNoBath : public Polarizability {
public:
  // constructors
  PolarizabilityNoBath(double omega_a, double alpha_zero,
                       GreensTensor *greens_tensor);
  PolarizabilityNoBath(std::string input_file);

  // calculate the polarizability tensor
  void calculate_tensor(cx_mat::fixed<3, 3> &alpha,
                        Options_Polarizability opts);
};

#endif // POLARIZABILITYNOBATH_H
