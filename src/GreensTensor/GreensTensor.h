#ifndef GREENSTENSOR_H
#define GREENSTENSOR_H

#include <armadillo>
#include "../Calculations/Integrations.h"
using namespace arma;

struct Options;
// A Greens tensor class

/*!
* This is an abstract class that implements an isotropic and reciprocal Greens tensor.
*/
class GreensTensor
{
protected:
  double v, za, beta;

public:
  virtual void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega) =0;
  virtual void integrate_k_2d(cx_mat::fixed<3,3>& GT, double omega, double k, Options opts)  =0;
  virtual void integrate_k_1d(cx_mat::fixed<3,3>& GT, double omega, Options opts) =0;
  double get_v(){return this->v;}
  double get_za(){return this->za;}
  double get_beta(){return this->beta;}

};

// A struct for integration options struct 
struct Options
{
  bool fancy_R = false;
  bool fancy_I = false;
  bool fancy_I_kv = false;
  bool fancy_I_temp = false;
  bool fancy_I_kv_temp = false;
  arma::vec::fixed<2> indices = {-1,-1};
  double omega = NAN;
  double k = NAN;
  GreensTensor* class_pt;
};

#endif // GREENSTENSOR_H
