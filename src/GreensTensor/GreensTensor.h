#ifndef GREENSTENSOR_H
#define GREENSTENSOR_H

#include <string>
#include <armadillo>
#include "../Calculations/Integrations.h"
using namespace arma;

//A struct with integration options
struct Options_GreensTensor;
// A Greens tensor class

/*!
* This is an abstract class that implements an isotropic and reciprocal Greens tensor.
*/
class GreensTensor
{
protected:
  double v, za, beta;

public:

  // constructor
  GreensTensor(double v, double za, double beta):v(v), za(za), beta(beta) {};
  GreensTensor(std::string input_file);

  // calculate the whole Green's tensor
  virtual void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega) =0;

  // integrate over a two-dimensional k space
  virtual void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)  =0;

  // integrate over a one-dimensional k space
  virtual void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts) =0;

  // getter functions
  double get_v() const {return this->v;}
  double get_za() const {return this->za;}
  double get_beta() const {return this->beta;}

};

// A struct for integration options
struct Options_GreensTensor
{
 // Different options for the integrand
  bool fancy_R = false;
  bool fancy_I = false;
  bool fancy_I_kv = false;
  bool fancy_I_temp = false;
  bool fancy_I_kv_temp = false;

  // Indices of the 3x3 GreensTensor
  vec::fixed<2> indices = {-1,-1};

  // Value of omega for the integration of the k-Variables
  double omega = NAN;

  // k-vector for the omega integration
  vec::fixed<2>  kvec = {NAN,NAN};

  // Pointer to the GreensTensor to be able to access the attributes of the class eventhough the integrand is static
  GreensTensor* class_pt;
};

#endif // GREENSTENSOR_H
