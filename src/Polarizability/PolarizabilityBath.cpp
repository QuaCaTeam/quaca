#include "PolarizabilityBath.h"

PolarizabilityBath::PolarizabilityBath(double omega_a, double alpha_zero,
                                       MemoryKernel *mu,
                                       GreensTensor *greens_tensor)
    : Polarizability(omega_a, alpha_zero, greens_tensor), mu(mu) {}

PolarizabilityBath::PolarizabilityBath(const std::string &input_file)
    : Polarizability(input_file) {
  this->mu =
      MemoryKernelFactory::create(input_file, "Polarizability.MemoryKernel");
}

void PolarizabilityBath::calculate_tensor(cx_mat::fixed<3, 3> &alpha,
                                          Options_Polarizability opts) {
  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  double omega = opts.omega;

  // calculate diagonal entries
  cx_mat::fixed<3, 3> diag(fill::zeros);
  diag(0, 0) = diag(1, 1) = diag(2, 2) =
      omega_a * omega_a - omega * omega - I * omega * mu->mu(omega);

  // calculate integral over green's tensor with fancy R
  cx_mat::fixed<3, 3> greens_R;
  struct Options_GreensTensor opts_R;
  opts_R.fancy_complex = RE;
  opts_R.omega = omega;
  opts_R.class_pt = this->greens_tensor;
  this->greens_tensor->integrate_k(greens_R, opts_R);

  // calculate integral over green's tensor with fancy I
  cx_mat::fixed<3, 3> greens_I;
  struct Options_GreensTensor opts_I;
  opts_I.fancy_complex = IM;
  opts_I.omega = omega;
  opts_I.class_pt = this->greens_tensor;

  this->greens_tensor->integrate_k(greens_I, opts_I);

  // put everything together
  alpha =
      alpha_zero * omega_a * omega_a *
      inv(diag - alpha_zero * omega_a * omega_a * (greens_R + I * greens_I));

  if (opts.fancy_complex == IM) {
    alpha = (alpha - trans(alpha)) /
            (2.0 * I); // trans is hermitean conjugation in armadillo
  } else if (opts.fancy_complex == RE) {
    alpha = (alpha + trans(alpha)) /
            (2.0); // trans is hermitean conjugation in armadillo
  }
}
