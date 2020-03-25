#ifndef GREENSTENSORPLATEMAGNETIC_H
#define GREENSTENSORPLATEMAGNETIC_H

#include "GreensTensorPlate.h"
#include <armadillo>
#include <assert.h>

struct Options_GreensTensorMagnetic : Options_GreensTensor {
  Tensor_Options EB = IGNORE;
  Tensor_Options BE = IGNORE;
  Tensor_Options BB = IGNORE;
};

class GreensTensorPlateMagnetic : public GreensTensorPlate {
public:
  // constructors
  GreensTensorPlateMagnetic(std::string input_file);
  GreensTensorPlateMagnetic(double v, double za, double beta,
                            ReflectionCoefficients *reflection_coefficients,
                            double delta_cut, vec::fixed<2> rel_err);

  void calculate_tensor(cx_mat::fixed<3, 3> &GT,
                        Options_GreensTensor opts);

  void integrate_k(cx_mat::fixed<3, 3> &GT,
                                            Options_GreensTensorMagnetic opts); 
  void integrate_k(cx_mat::fixed<3, 3> &GT,
                                            Options_GreensTensor opts); 

  // integrands
  static double integrand_2d_k_magnetic_R(double kappa_double, void *opts);
  static double integrand_2d_k_magnetic_I(double kappa_double, void *opts);
  static double integrand_1d_k_magnetic_R(double phi, void *opts);
  static double integrand_1d_k_magnetic_I(double phi, void *opts);
};

#endif
