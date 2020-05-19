#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include "GreensTensor.h"
#include <cmath>
#include <complex>

class GreensTensorVacuum : public GreensTensor {
private:
  double relerr; //integration error along the k_v direction

public:
  // constructors
  GreensTensorVacuum(double v, double beta, double relerr);
  GreensTensorVacuum(const std::string& input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(double omega, vec::fixed<2> k,
                        cx_mat::fixed<3, 3> &GT) const override ;

  // integrate over a two-dimensional k space
  void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                   Tensor_Options fancy_complex,
                   Weight_Options weight_function) const override;

  // integrand for integration over one-dimensional k space
  double integrand_k(double kv, double omega, const vec::fixed<2> &indices,
                                         Tensor_Options fancy_complex,
                                         Weight_Options weight_function) const;

  double omega_ch() const override;
  double get_relerr() const { return this->relerr; };

  // print info
  void print_info(std::ostream &stream) const override;
};

#endif // GREENSTENSORVACUUM_H
