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
  GreensTensorVacuum(const std::string& input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a two-dimensional k space
  void integrate_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrand for integration over one-dimensional k space
  static double integrand_k(double k, void *opts);
  double omega_ch();
  double get_relerr() const { return this->relerr; };

  void print_info(std::ofstream &file);
};

#endif // GREENSTENSORVACUUM_H
