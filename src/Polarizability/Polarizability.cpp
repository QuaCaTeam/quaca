// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

#include "../GreensTensor/GreensTensorFactory.h"
#include "Polarizability.h"

Polarizability::Polarizability(double omega_a, double alpha_zero,
                               std::shared_ptr<GreensTensor> greens_tensor)
    : omega_a(omega_a), alpha_zero(alpha_zero),
      greens_tensor(std::move(greens_tensor)) {}

Polarizability::Polarizability(const std::string &input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // read parameters
  this->omega_a = root.get<double>("Polarizability.omega_a");
  this->alpha_zero = root.get<double>("Polarizability.alpha_zero");

  // read greens tensor
  this->greens_tensor = GreensTensorFactory::create(input_file);
}

double Polarizability::integrand_omega(double omega,
                                       const vec::fixed<2> &indices,
                                       Tensor_Options fancy_complex) const {
  cx_mat::fixed<3, 3> alpha;
  calculate_tensor(omega, alpha, fancy_complex);

  assert(imag(alpha(indices(0), indices(1))) == 0);
  double result = real(alpha(indices(0), indices(1)));

  return result;
}

double Polarizability::integrate_omega(const vec::fixed<2> &indices, Tensor_Options fancy_complex, double omega_min,
                                       double omega_max, double relerr,
                                       double abserr) const {
  auto F = [=](double x) -> double {
    return this->integrand_omega(x, indices, fancy_complex);
  };
  return cquad(F, omega_min, omega_max, relerr, abserr);
}
