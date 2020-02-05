#include "PowerSpectrumHarmOsc.h"
#include <armadillo>
#include <string>

using namespace arma;

// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(std::string input_file)
    : PowerSpectrum(input_file) {

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("PowerSpectrum.type");
  assert(type == "harmonic oscillator");
};

PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(GreensTensor *greens_tensor,
                                           Polarizability *polarizability)
    : PowerSpectrum(greens_tensor, polarizability){};

void PowerSpectrumHarmOsc::calculate(cx_mat::fixed<3, 3> &powerspectrum,
                                     Options_PowerSpectrum opts) {
  double omega = opts.omega;
  if (opts.full_spectrum) {
    cx_mat::fixed<3, 3> green(fill::zeros);
    Options_GreensTensor opts_g;
    opts_g.fancy_I_temp = true;
    opts_g.omega = omega;
    opts_g.class_pt = this->greens_tensor;
    // Compute the Green's tensor
    this->greens_tensor->integrate_1d_k(green, opts_g);

    cx_mat::fixed<3, 3> alpha(fill::zeros);
    Options_Polarizability opts_alpha;
    opts_alpha.omega = omega;
    // Compute the polarizability
    this->polarizability->calculate_tensor(alpha, opts_alpha);

    // Combine the previous results to a power spectrum, see eq. [3.9] in
    // Marty's PhD thesis
    powerspectrum = 1. / M_PI * alpha * green * trans(alpha);
  }
  if (opts.non_LTE) {
    cx_mat::fixed<3, 3> green(fill::zeros);
    Options_GreensTensor opts_g;
    opts_g.fancy_I_non_LTE = true;
    opts_g.omega = omega;
    opts_g.class_pt = this->greens_tensor;
    // Compute the Green's tensor
    this->greens_tensor->integrate_1d_k(green, opts_g);

    cx_mat::fixed<3, 3> alpha(fill::zeros);
    Options_Polarizability opts_alpha;
    opts_alpha.omega = omega;

    // Compute the polarizability
    this->polarizability->calculate_tensor(alpha, opts_alpha);

    // Combine the previous results to a power spectrum, see eq. [3.9] in
    // Marty's PhD thesis
    powerspectrum = alpha * green * trans(alpha);
  }
};
