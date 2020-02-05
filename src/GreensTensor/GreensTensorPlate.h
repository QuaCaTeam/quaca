#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include "GreensTensor.h"
#include "Permittivity/PermittivityFactory.h"
#include "ReflectionCoefficients/ReflectionCoefficientsFactory.h"
#include <armadillo>
#include <assert.h>
#include <cmath>
#include <complex>
// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

//! The class of the Green's tensor above a flat macroscopic surface
class GreensTensorPlate : public GreensTensor {
private:
  // reflection coefficients are needed to describe the surface's response
  ReflectionCoefficients *reflection_coefficients;
  // kappa_cut defines the numerical cut-off of the kappa integration
  double delta_cut;
  vec::fixed<2> rel_err = {NAN, NAN};
  double za;

public:
  // constructor with direct input
  GreensTensorPlate(double v, double za, double beta,
                    ReflectionCoefficients *reflection_coefficients, double delta_cut,
                    vec::fixed<2> rel_err)
      : GreensTensor(v, beta) {
    this->za = za;
    this->delta_cut = delta_cut;
    this->rel_err = rel_err;
    this->reflection_coefficients = reflection_coefficients;
  };

  // constructor with .ini file
  GreensTensorPlate(std::string input_file) : GreensTensor(input_file) {
    this->reflection_coefficients = ReflectionCoefficientsFactory::create(input_file);

    // Create a root
    pt::ptree root;

    // Load the ini file in this ptree
    pt::read_ini(input_file, root);

    // read parameters
    this->za = root.get<double>("GreensTensor.za");
    this->delta_cut = root.get<double>("GreensTensor.delta_cut");
    this->rel_err(0) = root.get<double>("GreensTensor.rel_err_0");
    this->rel_err(1) = root.get<double>("GreensTensor.rel_err_1");
  };
  // calculate the pure Green's tensor of this class
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);
  // integrates the Green's tensor with respect to kappa
  void integrate_2d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);
  // integrates the Green's tensor with respect to phi and kappa
  void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);
  // integrand of Green's tensor element with respect to the phi integration
  static double integrand_1d_k(double kx, void *opts);
  // integrand of Green's tensor element with respect to the kappa integration
  static double integrand_2d_k(double ky, void *opts);
  // getter functions
  std::complex<double> get_r_p(double omega, double k);
  std::complex<double> get_r_s(double omega, double k);
  double get_za() { return this->za; }
  double get_delta_cut() { return this->delta_cut; }
  double get_rel_err_0() { return this->rel_err(0); }
  double get_rel_err_1() { return this->rel_err(1); }
};

#endif // GREENSTENSORPLATE_H
