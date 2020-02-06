#include "QuantumFriction.h"
#include "../GreensTensor/GreensTensorFactory.h"
#include "../Polarizability/PolarizabilityFactory.h"
#include "../PowerSpectrum/PowerSpectrumFactory.h"
#include <armadillo>
#include <algorithm>
#include <vector>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace arma;

QuantumFriction::QuantumFriction(std::string input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);
  this-> relerr_omega = root.get<double>("Friction.relerr_omega");
  // read greens tensor
  this->greens_tensor = GreensTensorFactory::create(input_file);
  this->polarizability = PolarizabilityFactory::create(input_file);
  this->powerspectrum = PowerSpectrumFactory::create(input_file);
};

QuantumFriction::QuantumFriction(GreensTensor *greens_tensor,
                                 Polarizability *polarizability,
                                 PowerSpectrum *powerspectrum, double relerr_omega)
    : greens_tensor(greens_tensor), polarizability(polarizability),
      powerspectrum(powerspectrum), relerr_omega(relerr_omega){};

double QuantumFriction::calculate(Options_Friction opts) {
  double result;
  double omega_a =this->polarizability->get_omega_a();
  // Collect all specifically relevant point within the integration
  std::vector<double> lim = {0., 0.99*omega_a , 1.001*omega_a , this->greens_tensor->omega_ch()};
  // Sort the points and erase duplicates
  std::sort(lim.begin(), lim.end());
  auto last = std::unique(lim.begin(), lim.end()); 
  lim.erase(last, lim.end());

  // Start integration
  result = 0.;

  for (int i =0; i < lim.size()-1 ; i++){ 
  result += cquad(&friction_integrand, &opts, lim[i], lim[i+1], relerr_omega, std::abs(result)*relerr_omega);
  
  };
  // Perform last integration from the last significant point to infinity
  result += qagiu(&friction_integrand, &opts, lim[lim.size()-1], relerr_omega, std::abs(result) * relerr_omega);
  return result;
};

double QuantumFriction::friction_integrand(double omega, void *opts) {
  // Implementation of the integrand of eq. (4.3) in Marty's PhD thesis
  // Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  Options_Friction *opts_pt = static_cast<Options_Friction *>(opts);
  if (opts_pt->full_spectrum) {
    cx_mat::fixed<3, 3> green_kv(fill::zeros);
    cx_mat::fixed<3, 3> green_temp_kv(fill::zeros);
    cx_mat::fixed<3, 3> alpha_I(fill::zeros);
    cx_mat::fixed<3, 3> powerspectrum(fill::zeros);

    // computation of the Green's tensor in the first term of eq. (4.3)
    Options_GreensTensor opts_g;
    opts_g.omega = omega;
    opts_g.class_pt = opts_pt->class_pt->greens_tensor;

    opts_g.fancy_I_kv = true;
    opts_pt->class_pt->greens_tensor->integrate_1d_k(green_kv, opts_g);

    // computation of the Green's tensor in the second term of eq. (4.3)
    opts_g.fancy_I_kv = false;
    opts_g.fancy_I_kv_temp = true;
    opts_pt->class_pt->greens_tensor->integrate_1d_k(green_temp_kv, opts_g);

    // computation of the imaginary part of the polarizability appearing in the
    // second term of eq. (4.3)
    Options_Polarizability opts_alpha;
    opts_alpha.omega = omega;
    opts_alpha.fancy_I = true;
    opts_pt->class_pt->polarizability->calculate_tensor(alpha_I, opts_alpha);

    // computation of the powerspectrum apperaing in the first term of eq. (4.3)
    Options_PowerSpectrum opts_S;
    opts_S.omega = omega;
    opts_S.full_spectrum = true;
    opts_pt->class_pt->powerspectrum->calculate(powerspectrum, opts_S);
    return real(2. * trace(-powerspectrum * green_kv +
                           1. / M_PI * alpha_I * green_temp_kv));
    // return real(2.*trace(-powerspectrum*green_kv));
    //  return real(trace(powerspectrum)*pow(omega,2));
    // return real(trace( 1./M_PI*alpha*green_temp_kv));
  } else if (opts_pt->non_LTE) {
    cx_mat::fixed<3, 3> J(fill::zeros);
    cx_mat::fixed<3, 3> green_kv(fill::zeros);
    cx_mat::fixed<3, 3> alpha_fancy_I(fill::zeros);
    cx_mat::fixed<3, 3> green_fancy_I_kv_non_LTE(fill::zeros);

    Options_GreensTensor opts_g;
    opts_g.omega = omega;
    opts_g.class_pt = opts_pt->class_pt->greens_tensor;

    opts_g.fancy_I_kv = true;
    opts_pt->class_pt->greens_tensor->integrate_1d_k(green_kv, opts_g);

    opts_g.fancy_I_kv = false;
    opts_g.fancy_I_kv_non_LTE = true;
    opts_pt->class_pt->greens_tensor->integrate_1d_k(green_fancy_I_kv_non_LTE,
                                                     opts_g);

    Options_PowerSpectrum opts_J;
    opts_J.non_LTE = true;
    opts_J.omega = omega;
    opts_pt->class_pt->powerspectrum->calculate(J, opts_J);

    Options_Polarizability opts_alpha;
    opts_alpha.fancy_I = true;
    opts_alpha.omega = omega;
    opts_pt->class_pt->polarizability->calculate_tensor(alpha_fancy_I,
                                                        opts_alpha);
    return real(
        trace(2. / M_PI *
              (-J * green_kv + alpha_fancy_I * green_fancy_I_kv_non_LTE)));
  } else {
    std::cerr << "No valid option for the calculation of quantum friction have "
                 "been passed"
              << std::endl;
    exit(0);
  }
  return 0;
};
