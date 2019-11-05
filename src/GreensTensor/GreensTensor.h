#ifndef GREENSTENSOR_H
#define GREENSTENSOR_H

#include <armadillo>
using namespace arma;

// A Greens tensor class
/*!
* This is an abstract class that implements an isotropic and reciprocal Greens tensor.
*/
class GreensTensor
{
protected:
  cx_mat::fixed<3,3> greens;
  double v, za;

public:
  virtual cx_mat::fixed<3,3> calculate_pure(double k_v, double omega) =0;
  virtual cx_mat::fixed<3,3> calculate_integrated(int flags, double omega) =0;

};

#endif // GREENSTENSOR_H
