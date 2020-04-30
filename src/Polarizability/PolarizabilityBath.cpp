#include "PolarizabilityBath.h"

#include <utility>

PolarizabilityBath::PolarizabilityBath(
    double omega_a, double alpha_zero, std::shared_ptr<MemoryKernel> mu,
    std::shared_ptr<GreensTensor> greens_tensor)
    : Polarizability(omega_a, alpha_zero, std::move(greens_tensor)),
      mu(std::move(mu)) {}

PolarizabilityBath::PolarizabilityBath(const std::string &input_file)
    : Polarizability(input_file) {
  this->mu =
      MemoryKernelFactory::create(input_file, "Polarizability.MemoryKernel");
}

void PolarizabilityBath::calculate_tensor(double omega,
                                          cx_mat::fixed<3, 3> &alpha,
                                          Tensor_Options fancy_complex) const {
  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // calculate diagonal entries
  cx_mat::fixed<3, 3> diag;
  diag.zeros();
  diag(0, 0) = diag(1, 1) = diag(2, 2) =
      omega_a * omega_a - omega * omega - I * omega * mu->calculate(omega);

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
