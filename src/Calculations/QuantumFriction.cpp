#include "QuantumFriction.h"
#include "../GreensTensor/GreensTensorFactory.h"
#include "../Polarizability/PolarizabilityFactory.h"
#include "../PowerSpectrum/PowerSpectrumFactory.h"
#include <armadillo>

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace arma;

QuantumFriction::QuantumFriction(std::string input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read greens tensor
  this->greens_tensor = GreensTensorFactory::create(input_file);
  this->polarizability = PolarizabilityFactory::create(input_file);
  this->powerspectrum = PowerSpectrumFactory::create(input_file);
};

QuantumFriction::QuantumFriction(GreensTensor *greens_tensor,
                                 Polarizability *polarizability,
                                 PowerSpectrum *powerspectrum)
    : greens_tensor(greens_tensor), polarizability(polarizability),
      powerspectrum(powerspectrum){};

double QuantumFriction::calculate(Options_Friction opts, double omega_min,
                                  double omega_max, double relerr,
                                  double epsabs) {
  return cquad(&friction_integrand, &opts, omega_min, omega_max, relerr,
               epsabs);
};

double QuantumFriction::friction_integrand(double omega, void *opts) {
  // Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  Options_Friction *opts_pt = static_cast<Options_Friction *>(opts);
  cx_mat::fixed<3, 3> green_kv(fill::zeros);
  cx_mat::fixed<3, 3> green_temp_kv(fill::zeros);
  cx_mat::fixed<3, 3> alpha(fill::zeros);
  cx_mat::fixed<3, 3> powerspectrum(fill::zeros);

  Options_GreensTensor opts_g;
  opts_g.omega = omega;
  opts_g.class_pt = opts_pt->class_pt->greens_tensor;

  opts_g.fancy_I_kv = true;
  opts_pt->class_pt->greens_tensor->integrate_1d_k(green_kv, opts_g);

  opts_g.fancy_I_kv = false;
  opts_g.fancy_I_kv_temp = true;
  opts_pt->class_pt->greens_tensor->integrate_1d_k(green_temp_kv, opts_g);

  // std::cout << "Greens tensor temp" << std::endl << green_temp_kv;
  // std::cout << "Greens tensor kv" << std::endl << green_kv;

  Options_Polarizability opts_alpha;
  opts_alpha.omega = omega;
  opts_pt->class_pt->polarizability->calculate_tensor(alpha, opts_alpha);
  // std::cout << "Polarizability" << std::endl << alpha;

  opts_pt->class_pt->powerspectrum->calculate(powerspectrum, omega);
  //  std::cout << "Power spectrum" << std::endl << powerspectrum;

  return real(2. * trace(-powerspectrum * green_kv + alpha * green_temp_kv));
};
