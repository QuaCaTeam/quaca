#ifndef GREENSTENSOR_H
#define GREENSTENSOR_H

#include "../Calculations/Integrations.h"
#include <armadillo>
using namespace arma;

// A struct with integration options
struct Options_GreensTensor;
// A Greens tensor class

/*!
 * This is an abstract class that implements an isotropic and reciprocal Greens
 * tensor.
 */
class GreensTensor {
protected:
  double v;
  double beta;

public:
  // constructor
  GreensTensor(double v, double beta) : v(v), beta(beta){};
  GreensTensor(std::string input_file);

  virtual void calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                Options_GreensTensor opts) = 0;
  virtual void integrate_2d_k(cx_mat::fixed<3, 3> &GT,
                              Options_GreensTensor opts) = 0;
  virtual void integrate_1d_k(cx_mat::fixed<3, 3> &GT,
                              Options_GreensTensor opts) = 0;
  // getter functions
  double get_v() { return this->v; }
  double get_beta() { return this->beta; }
  // calculate characteristic frequencies
  virtual double omega_ch() = 0;
};

// A struct for integration options
struct Options_GreensTensor {
  // Different options for the integrand
  bool fancy_R = false;
  bool fancy_I = false;
  bool fancy_I_kv = false;
  bool fancy_I_temp = false;
  bool fancy_I_non_LTE = false;
  bool fancy_I_kv_temp = false;
  bool fancy_I_kv_non_LTE = false;
  // Indices of the 3x3 GreensTensor
  arma::vec::fixed<2> indices = {-1, -1};
  // Value of omega for the integration of the k-Variables
  double omega = NAN;
  // k-vector for the omega integration
  vec::fixed<2> kvec = {NAN, NAN};
  // Pointer to the GreensTensor to be able to access the attributes of the
  // class eventhough the integrand is static
  GreensTensor *class_pt;
};

#endif // GREENSTENSOR_H
