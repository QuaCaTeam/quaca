#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include <complex>
#include <cmath>
#include "GreensTensor.h"

class GreensTensorVacuum : public GreensTensor
{

public:

  // constructors
  GreensTensorVacuum(double v, double beta): GreensTensor(v, 0, beta) {};
  GreensTensorVacuum(std::string input_file): GreensTensor(input_file) {};

  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  static double integrand_1d_k(double k, void* opts);

};


#endif // GREENSTENSORVACUUM_H
