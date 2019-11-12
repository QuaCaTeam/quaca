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
  double v, za;

public:
  virtual void calculate_pure(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega) =0;
  virtual void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, double kv, std::vector<std::string> options) =0;
  virtual void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, std::vector<std::string> options) =0;

};

#endif // GREENSTENSOR_H
