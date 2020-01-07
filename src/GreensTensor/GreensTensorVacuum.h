#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include <complex>
#include <cmath>
#include "GreensTensor.h"

class GreensTensorVacuum : public GreensTensor
{
  private:
    double za;

public:

  GreensTensorVacuum(double v, double beta);
  void calculate_tensor(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  static double integrand_1d_k(double k, void* opts);
  double get_za(){return this->za;};
  void set_za(double za){this->za = za;};

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  // integrand for integration over one-dimensional k space
  static double integrand_1d_k(double k, void* opts);
};


#endif // GREENSTENSORVACUUM_H
