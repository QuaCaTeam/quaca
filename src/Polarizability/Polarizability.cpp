// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

#include "../GreensTensor/GreensTensorFactory.h"
#include "../MemoryKernel/MemoryKernelFactory.h"
#include "Polarizability.h"

// Constructor without an internal bath
Polarizability::Polarizability(double omega_a, double alpha_zero,
                               std::shared_ptr<GreensTensor> greens_tensor)
    : omega_a(omega_a), alpha_zero(alpha_zero),
      greens_tensor(std::move(greens_tensor)) {
  this->mu = nullptr;
}
// Constructor with internal bath mu
Polarizability::Polarizability(double omega_a, double alpha_zero,
                               std::shared_ptr<MemoryKernel> mu,
                               std::shared_ptr<GreensTensor> greens_tensor)
    : omega_a(omega_a), alpha_zero(alpha_zero), mu(std::move(mu)),
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

  // Check if the category Polarizability.MemoryKernel exists
  boost::optional<pt::ptree &> child =
      root.get_child_optional("Polarizability.MemoryKernel");

  if (child) {
    this->mu =
        MemoryKernelFactory::create(input_file, "Polarizability.MemoryKernel");
  } else {
    this->mu = nullptr;
  };
}

void Polarizability::calculate_tensor(double omega, cx_mat::fixed<3, 3> &alpha,
                                      Tensor_Options fancy_complex) const {
  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // calculate diagonal entries
  cx_mat::fixed<3, 3> diag;
  diag.zeros();
  if (mu != nullptr) {
    diag(0, 0) = diag(1, 1) = diag(2, 2) =
        omega_a * omega_a - omega * omega - I * omega * mu->calculate(omega);
  } else {
    diag(0, 0) = diag(1, 1) = diag(2, 2) = omega_a * omega_a - omega * omega;
  };
  // calculate integral over green's tensor with fancy R
  cx_mat::fixed<3, 3> greens_R;
  this->greens_tensor->integrate_k(omega, greens_R, RE, UNIT);

  // calculate integral over green's tensor with fancy I
  cx_mat::fixed<3, 3> greens_I;
  this->greens_tensor->integrate_k(omega, greens_I, IM, UNIT);

  // put everything together
  alpha =
      alpha_zero * omega_a * omega_a *
      inv(diag - alpha_zero * omega_a * omega_a * (greens_R + I * greens_I));

  if (fancy_complex == IM) {
    alpha = (alpha - trans(alpha)) /
            (2.0 * I); // trans is hermitean conjugation in armadillo
  } else if (fancy_complex == RE) {
    alpha = (alpha + trans(alpha)) /
            (2.0); // trans is hermitean conjugation in armadillo
  }
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

double Polarizability::integrate_omega(const vec::fixed<2> &indices,
                                       Tensor_Options fancy_complex,
                                       double omega_min, double omega_max,
                                       double relerr, double abserr) const {
  auto F = [=](double x) -> double {
    return this->integrand_omega(x, indices, fancy_complex);
  };
  return cquad(F, omega_min, omega_max, relerr, abserr);
}
