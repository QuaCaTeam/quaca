#include "PowerSpectrumHarmOsc.h"
#include "../Polarizability/PolarizabilityBath.h"
#include <armadillo>
#include <cassert>
#include <string>

using namespace arma;

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

// Constructor with json-file
PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(const std::string& input_file)
    : PowerSpectrum(input_file) {

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("PowerSpectrum.type");
  // read the type of the polarizability and set the bool to indicate,
  // whether the polarizablity has an interal bath
  std::string polarizability_type =
      root.get<std::string>("Polarizability.type");
  has_bath = polarizability_type == "bath";


  // Ensure that correct type has been chosen
  assert(type == "harmonic oscillator");
}

// Constructor with initialization list
PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(
    std::shared_ptr<GreensTensor> greens_tensor,
    std::shared_ptr<Polarizability> polarizability)
    : PowerSpectrum(std::move(greens_tensor), std::move(polarizability)) {

  // read the type of the polarizability and set the bool to indicate,
  // whether the polarizablity has an interal bath
  auto pt = std::static_pointer_cast<PolarizabilityBath>(this->polarizability);
  has_bath = !(pt == nullptr);
}

// Compute the power spectrum for a given frequency \omega
void PowerSpectrumHarmOsc::calculate(cx_mat::fixed<3, 3> &powerspectrum,
                                     Options_PowerSpectrum opts) {
  // Read out variables double omega = opts.omega;
  double omega = opts.omega;

  // Compute the full spectrum
  if (opts.spectrum == FULL) {
    // Initialize tensor storing the Green's tensor and setting the integration
    // options for the Green's tensor
    cx_mat::fixed<3, 3> green(fill::zeros);
    Options_GreensTensor opts_g;
    opts_g.fancy_complex = IM;
    opts_g.weight_function = TEMP;
    opts_g.omega = omega;
    opts_g.class_pt = this->greens_tensor;

    // Compute the Green's tensor
    this->greens_tensor->integrate_k(green, opts_g);

    // Initialize tensor storing the polarizability and seting the integration
    // options for the polarizability
    cx_mat::fixed<3, 3> alpha(fill::zeros);
    Options_Polarizability opts_alpha;
    opts_alpha.omega = omega;
    if (opts_alpha.fancy_complex == IM)
      std::cout << "Is IM" << std::endl;
    if (opts_alpha.fancy_complex == RE)
      std::cout << "Is RE" << std::endl;

    // Compute the polarizability
    this->polarizability->calculate_tensor(alpha, opts_alpha);

    // Combine the Green's tensor and the polarizability, see eq. [3.9] in
    // Marty's PhD thesis
    powerspectrum = 1. / M_PI * alpha * green * trans(alpha);

    // check wether polarizability has an internal bath
    if (has_bath) {
      // To be able to use the attributes of PolarizabilityBath we have to
      // dynamically cast the pointer on the Polarizablity to a pointer on
      // PolarizablityBath
      auto pt =
          std::dynamic_pointer_cast<PolarizabilityBath>(this->polarizability);

      cx_mat::fixed<3, 3> bathTerm(fill::zeros);
      bathTerm = 1. / M_PI * alpha * 1. /
                 (pt->get_alpha_zero() * pow(pt->get_omega_a(), 2)) * omega *
                 real(pt->get_mu(omega)) /
                 (1. - exp(-this->greens_tensor->get_beta() * omega)) *
                 trans(alpha);
      powerspectrum += bathTerm;
    }
  }

  // Compute only the non-LTE contributions to the power spectrum
  if (opts.spectrum == NON_LTE_ONLY) {
    // Initialize tensor storing the Green's tensor and setting the integration
    // options for the Green's tensor
    cx_mat::fixed<3, 3> green(fill::zeros);
    Options_GreensTensor opts_g;
    opts_g.fancy_complex = IM;
    opts_g.weight_function = NON_LTE;
    opts_g.omega = omega;
    opts_g.class_pt = this->greens_tensor;

    // Compute the Green's tensor
    this->greens_tensor->integrate_k(green, opts_g);

    // Initialize tensor storing the polarizability and seting the integration
    // options for the polarizability
    cx_mat::fixed<3, 3> alpha(fill::zeros);
    Options_Polarizability opts_alpha;
    opts_alpha.omega = omega;

    // Compute the polarizability
    this->polarizability->calculate_tensor(alpha, opts_alpha);

    // Combine the Green's tensor and the polarizability, see eq. [3.9] in
    // Marty's PhD thesis
    powerspectrum = alpha * green * trans(alpha);
  }
}
