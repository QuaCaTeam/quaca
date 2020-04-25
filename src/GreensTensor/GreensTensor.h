#ifndef GREENSTENSOR_H
#define GREENSTENSOR_H

#include "../Calculations/Integrations.h"
#include <armadillo>
using namespace arma;

// Enum variables for the different integration options
enum Tensor_Options { COMPLEX, IM, RE };
enum Weight_Options { UNIT, KV, TEMP, KV_TEMP, NON_LTE, KV_NON_LTE };
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
  GreensTensor(double v, double beta);
  GreensTensor(const std::string& input_file);

  // calculate the tensor in frequency and momentum space
  virtual void calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                Options_GreensTensor opts) = 0;

  // integrate over a two-dimensional k space
  virtual void integrate_k(cx_mat::fixed<3, 3> &GT,
                           Options_GreensTensor opts) = 0;

  // getter functions
  double get_v() const { return this->v; };
  double get_beta() const { return this->beta; };

  virtual double omega_ch() = 0;

  // setter function
  virtual void set_v(double v_new) { this->v = v_new; };

  virtual void print_info(std::ofstream &file) = 0;
};

// A struct for integration options
struct Options_GreensTensor {
  // Different options for the integrand
  Tensor_Options fancy_complex = COMPLEX;
  Weight_Options weight_function = UNIT;
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
