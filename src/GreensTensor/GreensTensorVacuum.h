#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include <complex>
#include <cmath>
#include "GreensTensor.h"

class GreensTensorVacuum : public GreensTensor
{

public:

  GreensTensorVacuum(double v, double beta);
  GreensTensorVacuum(std::string input_file);
  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  static double integrand_k_1d(double k, void* opts);

};


#endif // GREENSTENSORVACUUM_H
