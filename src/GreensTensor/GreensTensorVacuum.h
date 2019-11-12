#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include <complex>
#include <cmath>
#include "GreensTensor.h"

class GreensTensorVacuum : public GreensTensor
{

public:

  GreensTensorVacuum(double v, double beta);
  void calculate_pure(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, double kv);
  void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, Options *opts);

};


#endif // GREENSTENSORVACUUM_H
