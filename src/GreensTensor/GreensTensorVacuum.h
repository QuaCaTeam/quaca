#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include "GreensTensor.h"
#include <cmath>
#include <complex>

class GreensTensorVacuum : public GreensTensor {
public:
  // constructors
  GreensTensorVacuum(double v, double beta);
  GreensTensorVacuum(std::string input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrand for integration over one-dimensional k space
  static double integrand_1d_k(double k, void *opts);
};

#endif // GREENSTENSORVACUUM_H
