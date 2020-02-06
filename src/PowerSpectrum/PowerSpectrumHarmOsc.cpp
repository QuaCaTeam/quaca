#include "PowerSpectrumHarmOsc.h"
#include <armadillo>
#include <string>

using namespace arma;

// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

//Constructor with ini-file
PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(std::string input_file)
    : PowerSpectrum(input_file) {

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("PowerSpectrum.type");

  //Ensure that correct type has been chosen
  assert(type == "harmonic oscillator");
};

//Constructor with initialization list
PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(GreensTensor *greens_tensor,
                                           Polarizability *polarizability)
    : PowerSpectrum(greens_tensor, polarizability){};

//Compute the power spectrum for a given frequency \omega
void PowerSpectrumHarmOsc::calculate(cx_mat::fixed<3, 3> &powerspectrum,
                                     Options_PowerSpectrum opts) {
    //Read out variables
  double omega = opts.omega;

  //Compute the full spectrum
  if (opts.full_spectrum) {

      //Initialize tensor storing the Green's tensor and setting the integration 
      //options for the Green's tensor
    cx_mat::fixed<3, 3> green(fill::zeros);
    Options_GreensTensor opts_g;
    opts_g.fancy_I_temp = true;
    opts_g.omega = omega;
    opts_g.class_pt = this->greens_tensor;

    // Compute the Green's tensor
    this->greens_tensor->integrate_1d_k(green, opts_g);

    //Initialize tensor storing the polarizability and seting the integration
    //options for the polarizability
    cx_mat::fixed<3, 3> alpha(fill::zeros);
    Options_Polarizability opts_alpha;
    opts_alpha.omega = omega;

    // Compute the polarizability
    this->polarizability->calculate_tensor(alpha, opts_alpha);

    // Combine the Green's tensor and the polarizability, see eq. [3.9] in
    // Marty's PhD thesis
    powerspectrum = 1. / M_PI * alpha * green * trans(alpha);
  }

  //Compute only the non-LTE contributions to the power spectrum
  if (opts.non_LTE) {

    //Initialize tensor storing the Green's tensor and setting the integration 
    //options for the Green's tensor
    cx_mat::fixed<3, 3> green(fill::zeros);
    Options_GreensTensor opts_g;
    opts_g.fancy_I_non_LTE = true;
    opts_g.omega = omega;
    opts_g.class_pt = this->greens_tensor;

    // Compute the Green's tensor
    this->greens_tensor->integrate_1d_k(green, opts_g);

    //Initialize tensor storing the polarizability and seting the integration
    //options for the polarizability
    cx_mat::fixed<3, 3> alpha(fill::zeros);
    Options_Polarizability opts_alpha;
    opts_alpha.omega = omega;

    // Compute the polarizability
    this->polarizability->calculate_tensor(alpha, opts_alpha);

    // Combine the Green's tensor and the polarizability, see eq. [3.9] in
    // Marty's PhD thesis
    powerspectrum = alpha * green * trans(alpha);
  }
};
