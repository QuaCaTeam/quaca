#include <algorithm>
#include <vector>

#include <armadillo>
using namespace arma;

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "../GreensTensor/GreensTensorFactory.h"
#include "Friction.h"

Friction::Friction(const std::string &input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);
  this->relerr_omega = root.get<double>("Friction.relerr_omega");
  this->sym_filter = root.get<std::string>("Friction.sym_filter");

  // read greens tensor
  this->powerspectrum = std::make_shared<PowerSpectrum>(input_file);
  this->polarizability = powerspectrum->get_polarizability();
  this->greens_tensor = powerspectrum->get_greens_tensor();
}

Friction::Friction(std::shared_ptr<GreensTensor> greens_tensor,
                   std::shared_ptr<Polarizability> polarizability,
                   std::shared_ptr<PowerSpectrum> powerspectrum,
                   double relerr_omega, std::string sym_filter)
    : greens_tensor(greens_tensor), polarizability(polarizability),
      powerspectrum(powerspectrum), relerr_omega(relerr_omega),
	sym_filter(sym_filter) {}

double Friction::calculate(Spectrum_Options spectrum) const {
  double result;
  double omega_a = this->polarizability->get_omega_a();
  // Collect all specifically relevant point within the integration
  std::vector<double> lim = {0., 0.99 * omega_a, 1.001 * omega_a,
                             this->greens_tensor->omega_ch()};
  // Sort the points and erase duplicates
  std::sort(lim.begin(), lim.end());
  auto last = std::unique(lim.begin(), lim.end());
  lim.erase(last, lim.end());

  // Start integration
  result = 0.;
  auto F = [=](double x) -> double { return friction_integrand(x, spectrum); };

  for (int i = 0; i < (int)lim.size() - 1; i++) {
    result += cquad(F, lim[i], lim[i + 1], relerr_omega,
                    std::abs(result) * relerr_omega);
  }
  // Perform last integration from the last significant point to infinity
  result += qagiu(F, lim[lim.size() - 1], relerr_omega,
                  std::abs(result) * relerr_omega);
  return result;
}

double Friction::friction_integrand(double omega,
                                    Spectrum_Options spectrum) const {
  // Compute the full spectrum of the power spectrum
  if (spectrum == FULL) {

    // Initialize all tensors
    cx_mat::fixed<3, 3> green_kv(fill::zeros);
    cx_mat::fixed<3, 3> green_temp_kv(fill::zeros);
    cx_mat::fixed<3, 3> alpha_I(fill::zeros);
    cx_mat::fixed<3, 3> powerspec(fill::zeros);

    // computation of the Green's tensor in the first term of eq. (4.3)
    greens_tensor->integrate_k(omega, green_kv, IM, KV);

    // computation of the Green's tensor in the second term of eq. (4.3)
    greens_tensor->integrate_k(omega, green_temp_kv, IM, KV_TEMP);

    // computation of the imaginary part of the polarizability appearing in the
    // second term of eq. (4.3)
    polarizability->calculate_tensor(omega, alpha_I, IM);

    // computation of the powerspectrum apperaing in the first term of eq. (4.3)
    powerspectrum->calculate(omega, powerspec, FULL);

    return real(2. * trace(-powerspec * green_kv +
                           1. / M_PI * alpha_I * green_temp_kv));
  }
  // Compute onyl the non-LTE contribution of the power-spectrum
  else if (spectrum == NON_LTE_ONLY) {
    // Initialize all tensors
    cx_mat::fixed<3, 3> J(fill::zeros);
    cx_mat::fixed<3, 3> green_kv(fill::zeros);
    cx_mat::fixed<3, 3> alpha_fancy_I(fill::zeros);
    cx_mat::fixed<3, 3> green_fancy_I_kv_non_LTE(fill::zeros);

    // Compute the Green's tensor in the first term of eq. (4.5)
    greens_tensor->integrate_k(omega, green_kv, IM, KV);

    // Compute the Green's tensor in the second term of eq. (4.5)
    // the \Sigma distribution is already included here
    greens_tensor->integrate_k(omega, green_fancy_I_kv_non_LTE, IM, KV_NON_LTE);

    // Compute the power spectrum for the first term of eq. (4.5)
    powerspectrum->calculate(omega, J, NON_LTE_ONLY);

    // Compute the polarizability for the second term of eq. (4.5)
    polarizability->calculate_tensor(omega, alpha_fancy_I, IM);

    // Optionally, apply filter with respect to matrix symmetries
    if (this->sym_filter == "symmetric") {
	    J(2, 0) = 0;
	    J(0, 2) = 0;
	    alpha_fancy_I(2, 0) = 0;
	    alpha_fancy_I(0, 2) = 0;
    }
    else if (this->sym_filter == "antisymmetric") {
	    J(0, 0) = 0;
	    J(1, 1) = 0;
	    J(2, 2) = 0;
	    alpha_fancy_I(0, 0) = 0;
	    alpha_fancy_I(1, 1) = 0;
	    alpha_fancy_I(2, 2) = 0;
    }

    // Put everything together
    return real(
        trace(2. *
              (-J * green_kv + alpha_fancy_I * green_fancy_I_kv_non_LTE)));
  }
  // Ensure that a valid integration argument has been passed
  else {
    std::cerr << "No valid option for the calculation of quantum friction have "
                 "been passed"
              << std::endl;
    exit(0);
  }
}

void Friction::print_info(std::ostream &stream) const {
  stream << "# Friction\n#\n"
  << "# relerr_omega = " << relerr_omega << "\n";
 greens_tensor->print_info(stream);
 polarizability->print_info(stream);
 powerspectrum->print_info(stream);
}
