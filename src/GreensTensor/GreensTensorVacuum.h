#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include "GreensTensor.h"
#include <cmath>
#include <complex>

class GreensTensorVacuum : public GreensTensor {
private:
  double relerr;

public:
  // constructors
  GreensTensorVacuum(double v, double beta, double relerr);
  GreensTensorVacuum(std::string input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(double omega, vec::fixed<2> k,
                        cx_mat::fixed<3, 3> &GT) const;

  // integrate over a two-dimensional k space
  void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                   Tensor_Options fancy_complex,
                   Weight_Options weight_function) const;

  // integrand for integration over one-dimensional k space
  double integrand_k(double kv, double omega, const vec::fixed<2> &indices,
                                         Tensor_Options fancy_complex,
                                         Weight_Options weight_function) const;

  double omega_ch() const;
  double get_relerr() const { return this->relerr; };
};

#endif // GREENSTENSORVACUUM_H
