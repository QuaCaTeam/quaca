#ifndef GREENSTENSORPLATEMAGNETIC_H
#define GREENSTENSORPLATEMAGNETIC_H

#include "GreensTensorPlate.h"
#include <armadillo>
#include <assert.h>

class GreensTensorPlateMagnetic : public GreensTensorPlate {
public:
  // constructors
  GreensTensorPlateMagnetic(std::string input_file);
  GreensTensorPlateMagnetic(double v, double beta, double za,
                            std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
                            double delta_cut, vec::fixed<2> rel_err);

  void calculate_tensor(double omega, vec::fixed<2> k, cx_mat::fixed<3, 3> &GT);

  void integrate_k(double omega, cx_mat::fixed<3, 3> &GT, Tensor_Options EE, Tensor_Options BE, Tensor_Options EB,
                   Tensor_Options BB, Weight_Options weight_function);

  // integrands
  double integrand_2d_k_R(double kappa_double, double omega, double phi, const vec::fixed<2> &indices,
                                 Tensor_Options EE, Tensor_Options EB, Tensor_Options BE, Tensor_Options BB,
                                 Weight_Options weight_function);
  double integrand_2d_k_I(double kappa_double, double omega, double phi, const vec::fixed<2> &indices,
                                 Tensor_Options EE, Tensor_Options EB, Tensor_Options BE, Tensor_Options BB,
                                 Weight_Options weight_function);
  double integrand_1d_k_R(double phi, double omega, const vec::fixed<2> &indices, Tensor_Options EE,
                                 Tensor_Options EB,
                                 Tensor_Options BE, Tensor_Options BB, Weight_Options weight_function);
  double integrand_1d_k_I(double phi, double omega, const vec::fixed<2> &indices, Tensor_Options EE,
                                 Tensor_Options EB,
                                 Tensor_Options BE, Tensor_Options BB, Weight_Options weight_function);
};

#endif
