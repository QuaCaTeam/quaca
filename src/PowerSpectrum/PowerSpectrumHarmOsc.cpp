#include "PowerSpectrumHarmOsc.h"
#include "../Polarizability/PolarizabilityBath.h"
#include <armadillo>
#include <cassert>
#include <string>

using namespace arma;

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

// Constructor with json-file
PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(const std::string &input_file)
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
  has_bath = (polarizability_type.compare("bath") == 0);

  // Ensure that correct type has been chosen
  assert(type == "harmonic oscillator");
}

// Constructor with initialization list
PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(
    std::shared_ptr<GreensTensor> greens_tensor,
    std::shared_ptr<Polarizability> polarizability)
    : PowerSpectrum(greens_tensor, polarizability) {

  // read the type of the polarizability and set the bool to indicate,
  // whether the polarizablity has an interal bath
  auto pt = std::dynamic_pointer_cast<PolarizabilityBath>(
      this->polarizability);
  has_bath = !(pt == nullptr);
}

// Compute the power spectrum for a given frequency \omega
void PowerSpectrumHarmOsc::calculate(double omega, cx_mat::fixed<3, 3> &powerspectrum,
                                     Spectrum_Options spectrum) const {
  // Compute the full spectrum
  if (spectrum == FULL) {
    // Initialize tensor storing the Green's tensor and setting the integration
    // options for the Green's tensor
    cx_mat::fixed<3, 3> green(fill::zeros);
    greens_tensor->integrate_k(omega, green, IM, TEMP);

    // Initialize tensor storing the polarizability and seting the integration
    // options for the polarizability
    cx_mat::fixed<3, 3> alpha(fill::zeros);
    polarizability->calculate_tensor(omega, alpha, COMPLEX);

    // Combine the Green's tensor and the polarizability, see eq. [3.9] in
    // Marty's PhD thesis
    powerspectrum = 1. / M_PI * alpha * green * trans(alpha);

    // check wether polarizability has an internal bath
    if (has_bath) {
      // To be able to use the attributes of PolarizabilityBath we have to
      // dynamically cast the pointer on the Polarizablity to a pointer on
      // PolarizablityBath
      auto pt = std::dynamic_pointer_cast<PolarizabilityBath>(
          polarizability);

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
  if (spectrum == NON_LTE_ONLY) {
    // Initialize tensor storing the Green's tensor and setting the integration
    // options for the Green's tensor
    cx_mat::fixed<3, 3> green(fill::zeros);

    // Compute the Green's tensor
    this->greens_tensor->integrate_k(omega, green, IM, NON_LTE);

    // Initialize tensor storing the polarizability and seting the integration
    // options for the polarizability
    cx_mat::fixed<3, 3> alpha(fill::zeros);
    // Compute the polarizability
    this->polarizability->calculate_tensor(omega, alpha, COMPLEX);

    // Combine the Green's tensor and the polarizability, see eq. [3.9] in
    // Marty's PhD thesis
    powerspectrum = alpha * green * trans(alpha);
  }
}
