#ifndef GREENSTENSOR_H
#define GREENSTENSOR_H

#include <armadillo>
using namespace arma;

// A struct for integration options
struct Options
{
  bool fancy_R = false;
  bool fancy_I = false;
  bool fancy_I_kv = false;
  bool fancy_I_temp = false;
  bool fancy_I_kv_temp = false;
};

// A Greens tensor class
/*!
* This is an abstract class that implements an isotropic and reciprocal Greens tensor.
*/
class GreensTensor
{
protected:
  double v, za, beta;

public:
  virtual void calculate_pure(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega) =0;
  virtual void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, Options *opts) =0;
  virtual void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, double kv) =0;

};

#endif // GREENSTENSOR_H
