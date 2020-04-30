#ifndef GREENSTENSOR_H
#define GREENSTENSOR_H

#include "../Calculations/Integrations.h"
#include <armadillo>
using namespace arma;

// Enum variables for the different integration options
enum Tensor_Options { COMPLEX, IM, RE };
enum Weight_Options { UNIT, KV, TEMP, KV_TEMP, NON_LTE, KV_NON_LTE };

//! A Greens tensor class
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
  explicit GreensTensor(const std::string &input_file);

  // calculate the tensor in frequency and momentum space
  virtual void calculate_tensor(double omega, vec::fixed<2> k,
                                cx_mat::fixed<3, 3> &GT) const = 0;

  // integrate over a two-dimensional k space
  virtual void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                           Tensor_Options fancy_complex,
                           Weight_Options weight_function) const = 0;

  // getter functions
  double get_v() const { return this->v; };
  double get_beta() const { return this->beta; };

  virtual double omega_ch() const = 0;

  // setter function
  virtual void set_v(double v_new) { this->v = v_new; };
};

#endif // GREENSTENSOR_H
