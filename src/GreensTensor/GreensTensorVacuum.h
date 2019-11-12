#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include <complex>
#include <cmath>
#include "GreensTensor.h"

class GreensTensorVacuum : public GreensTensor
{
protected:
  double v;

public:

  GreensTensorVacuum(double a);
  void calculate_pure(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, double kv, std::vector<std::string> options);
  void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, std::vector<std::string> options);

};


#endif // GREENSTENSORVACUUM_H
