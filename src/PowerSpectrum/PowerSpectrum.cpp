#include "PowerSpectrum.h"
#include "../GreensTensor/GreensTensorFactory.h"

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

// Constructor using a json-file
PowerSpectrum::PowerSpectrum(const std::string &input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // initialize polarizability by an input file
  this->polarizability = std::make_shared<Polarizability>(input_file);
  // initialize Green's tensor by an input file
  this->greens_tensor = polarizability->get_greens_tensor(); 
}

// Constructor with initialization list
PowerSpectrum::PowerSpectrum(std::shared_ptr<GreensTensor> greens_tensor,
                             std::shared_ptr<Polarizability> polarizability)
    : greens_tensor(std::move(greens_tensor)), polarizability(std::move(polarizability)) {}

// Compute the power spectrum for a given frequency \omega
void PowerSpectrum::calculate(double omega, cx_mat::fixed<3, 3> &powerspectrum,
                              Spectrum_Options spectrum) const {
  // imaginary unit
  std::complex<double> I(0.0, 1.0);
  // Compute the full spectrum
  if (spectrum == FULL) {
    // Initialize tensor storing the Green's tensor and setting the integration
    // options for the Green's tensor
    cx_mat::fixed<3, 3> green(fill::zeros);
    greens_tensor->integrate_k(omega, green, IM, NON_LTE);

    // Initialize tensor storing the polarizability and seting the integration
    // options for the polarizability
    cx_mat::fixed<3, 3> alpha(fill::zeros);
    polarizability->calculate_tensor(omega, alpha, COMPLEX);

    cx_mat::fixed<3, 3> alphaI(fill::zeros);
    alphaI = (alpha - trans(alpha)) / (2. * I);

    // First setting the LTE part
    powerspectrum =
        (alphaI / M_PI) / (1. - exp(-this->greens_tensor->get_beta() * omega));

    // Adding the NON-LTE part
    powerspectrum += 1. / M_PI * alpha * green * trans(alpha);
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
    powerspectrum = 1. / M_PI * alpha * green * trans(alpha);
  }
}

void PowerSpectrum::print_info(std::ostream &stream) const {
  stream << "# PowerSpectrum\n#\n";
  greens_tensor->print_info(stream);
  polarizability->print_info(stream);
}
