#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>
#include <iomanip> // std::setprecision

TEST_CASE("The tensors from calculate_tensor and integrand_2d_k coincide",
          "[GreensTensorPlateMagnetic]") {
  // Here we considered also the volume element from the integration.
  std::complex<double> I(0.0, 1.0);
  // Calculate computes all possible entries of the electric Green's tensor. On
  // the other hand integrand_2d_k_magnetic already sets all entries linear in
  // k_y to zero, as the integral from 0 to 2 Pi would vanish there phi has to
  // be set, such that sin(phi) = 0
  auto phi =  GENERATE(0., M_PI);
  auto k = GENERATE(take(3, random(0., 1e2)));
  auto omega = GENERATE(take(3, random(-1e2, 1e2)));
  double k_x = k * cos(phi);
  double k_y = k * sin(phi);
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  double kappa_double;
  std::complex<double> kappa, volume_element;
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic Greens(v, za, 0.1, &refl, delta_cut, rel_err);
  struct Options_GreensTensorMagnetic opts;
  opts.class_pt = &Greens;

  // First, the calculate_tensor operation is used to generate the
  // Green's tensor with fancy_I
  cx_mat::fixed<3, 3> Green(fill::zeros);
  cx_mat::fixed<3, 3> Green_fancy_I(fill::zeros);
  cx_mat::fixed<3, 3> Green_fancy_R(fill::zeros);

  opts.kvec(0) = k_x;
  opts.kvec(1) = k_y;
  opts.omega = omega + k_x * v;
  k = sqrt(k_x * k_x + k_y * k_y);
  kappa = sqrt(std::complex<double>(k * k - opts.omega * opts.omega, 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));
  volume_element = (real(kappa)-imag(kappa))*k / (k - cos(phi) * v * opts.omega);

  Greens.calculate_tensor(Green, opts);

  Green_fancy_I = (Green - trans(Green)) / (2. * I);
  Green_fancy_R = (Green + trans(Green)) / (2.);
  Green *= volume_element;
  Green_fancy_I *= volume_element;
  Green_fancy_R *= volume_element;

  //variable to store the kappa for the sigma matrix.
  std::complex<double> kappa_sigma = kappa;
  //ensure reality of kappa_sigma for complex kappa
  if(opts.omega < 0)
  {
    kappa_sigma = conj(kappa_sigma);

  }
  cx_mat::fixed<3, 3> Sigma = {
      {0, 0, 0}, {sin(phi), -cos(phi), 0}, {I * kappa_sigma / k, 0, -cos(phi)}};
  Sigma *= k * v / (opts.omega);

  cx_mat::fixed<3, 3> LHS(fill::zeros);
  cx_mat::fixed<3, 3> RHS(fill::zeros);

  if (kappa.real() == 0.) {
    kappa_double = kappa.imag();
  } else {
    kappa_double = kappa.real();
  }

  opts.kvec(0) = phi;
  opts.omega = omega;
  // Second, the integrand_2d_k operation is used
  SECTION("Options: EE, IM") {
    opts.fancy_complex = IM;
    // loop over all indices
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) =
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += (2 * M_PI) * I *
                     Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Green_fancy_I;
    REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
  }

  SECTION("Options: EE, RE") {
    opts.fancy_complex = RE;
    // loop over all indices
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) =
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += (2 * M_PI) * I *
                     Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Green_fancy_R;
    REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
  }
  opts.fancy_complex = IGNORE;

  SECTION("Options: BE, RE") {
    opts.BE = RE;
    // loop over all indices
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) =
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += (2 * M_PI) * I *
                     Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Sigma * Green_fancy_R;
    REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
  }

  SECTION("Options: BE, IM") {
    opts.BE = IM;
    // loop over all indices
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) =
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += (2 * M_PI) * I *
                     Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Sigma * Green_fancy_I;
    REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
  }
  opts.BE = IGNORE;

  SECTION("Options: EB, RE") {
    opts.EB = RE;
    // loop over all indices
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) =
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += (2 * M_PI) * I *
                     Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Green_fancy_R * trans(Sigma);
    REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
  }
  SECTION("Options: EB, IM") {
    opts.EB = IM;
    // loop over all indices
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) =
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += (2 * M_PI) * I *
                     Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Green_fancy_I * trans(Sigma);
    REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
  }

  opts.EB = IGNORE;

  SECTION("Options: BB, RE") {
    opts.BB = RE;
    // loop over all indices
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) =
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += (2 * M_PI) * I *
                     Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Sigma * Green_fancy_R * trans(Sigma);
    REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
  }
  SECTION("Options: BB, IM") {
    opts.BB = IM;
    // loop over all indices
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) =
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += (2 * M_PI) * I *
                     Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Sigma * Green_fancy_I * trans(Sigma);
    REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
  }
  opts.BB = IGNORE;
  

}


TEST_CASE("Crossing relation in k-space works for all Green's tensors",
          "[GreensTensorPlateMagnetic]") {
  
  // Here we considered also the volume element from the integration.
  std::complex<double> I(0.0, 1.0);
  // Calculate computes all possible entries of the electric Green's tensor. On
  // the other hand integrand_2d_k_magnetic already sets all entries linear in
  // k_y to zero, as the integral from 0 to 2 Pi would vanish there phi has to
  // be set, such that sin(phi) = 0
  auto phi = GENERATE(0., M_PI);
  auto k = GENERATE(take(3,random(0.,1e2)));
  auto omega = GENERATE(take(3,random(-1e2,1e2)));
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  double kappa_double;
  std::complex<double> kappa;
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic greens_tensor(v, za, 0.1, &refl, delta_cut,
                                          rel_err);
  struct Options_GreensTensorMagnetic opts;
  opts.class_pt = &greens_tensor;
  opts.omega = omega;
  opts.kvec(0) = phi;
  double omega_doppler = omega + k * cos(phi) * v;

  // Matrices to store the result, which we want to compare
  cx_mat::fixed<3, 3> LHS(fill::zeros);
  cx_mat::fixed<3, 3> RHS(fill::zeros);

  kappa = sqrt(std::complex<double>(k * k - pow(omega_doppler, 2), 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));

  if (kappa.real() == 0.) {
    kappa_double = kappa.imag();
  } else {
    kappa_double = kappa.real();
  }

  SECTION("G^EE_I fulfills the crossing relation in k-space") {
    opts.fancy_complex = IM;

    // Green's tensor with positive frequency and k-vector
    opts.omega = omega;
    opts.kvec(0) = phi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) += (2. * M_PI) *
                    greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += I * (2. * M_PI) *
                    greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    // Green's tensor with negative frequency and k-vector
    opts.omega = -omega;
    opts.kvec(0) = phi + M_PI;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        RHS(i, j) += (2. * M_PI) *
                    greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        RHS(i, j) += I * (2. * M_PI) *
                    greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    REQUIRE(approx_equal(conj(LHS), -RHS, "abs", 1e-12));
  }
  
  SECTION("G^EE_R fulfills the crossing relation in k-space") {
    opts.fancy_complex = RE;

    // Green's tensor with positive frequency and k-vector
    opts.omega = omega;
    opts.kvec(0) = phi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) += (2. * M_PI) *
                    greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    // Green's tensor with negative frequency and k-vector
    opts.omega = -omega;
    opts.kvec(0) = phi + M_PI;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        RHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        RHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    REQUIRE(approx_equal(conj(LHS), RHS, "abs", 1e-12));
  }
  opts.fancy_complex = IGNORE;
  SECTION("G^BE_I fulfills the crossing relation in k-space") {
    opts.BE = IM;

    // Green's tensor with positive frequency and k-vector
    opts.omega = omega;
    opts.kvec(0) = phi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    // Green's tensor with negative frequency and k-vector
    opts.omega = -omega;
    opts.kvec(0) = phi + M_PI;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        RHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        RHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    REQUIRE(approx_equal(conj(LHS), -RHS, "abs", 1e-12));
  }
  SECTION("G^BE_R fulfills the crossing relation in k-space") {
    opts.BE = RE;

    // Green's tensor with positive frequency and k-vector
    opts.omega = omega;
    opts.kvec(0) = phi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    // Green's tensor with negative frequency and k-vector
    opts.omega = -omega;
    opts.kvec(0) = phi + M_PI;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        RHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        RHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    REQUIRE(approx_equal(conj(LHS), RHS, "abs", 1e-12));
  }
  opts.BE = IGNORE;
  SECTION("G^EB_I fulfills the crossing relation in k-space") {
    opts.EB = IM;

    // Green's tensor with positive frequency and k-vector
    opts.omega = omega;
    opts.kvec(0) = phi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    // Green's tensor with negative frequency and k-vector
    opts.omega = -omega;
    opts.kvec(0) = phi + M_PI;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        RHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        RHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    REQUIRE(approx_equal(conj(LHS), -RHS, "abs", 1e-12));
  }
  SECTION("G^EB_R fulfills the crossing relation in k-space") {
    opts.EB = RE;

    // Green's tensor with positive frequency and k-vector
    opts.omega = omega;
    opts.kvec(0) = phi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    // Green's tensor with negative frequency and k-vector
    opts.omega = -omega;
    opts.kvec(0) = phi + M_PI;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        RHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        RHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    REQUIRE(approx_equal(conj(LHS), RHS, "abs", 1e-12));
  }
  opts.EB = IGNORE;
  SECTION("G^BB_I fulfills the crossing relation in k-space") {
    opts.BB = IM;

    // Green's tensor with positive frequency and k-vector
    opts.omega = omega;
    opts.kvec(0) = phi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    // Green's tensor with negative frequency and k-vector
    opts.omega = -omega;
    opts.kvec(0) = phi + M_PI;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        RHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        RHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    REQUIRE(approx_equal(conj(LHS), -RHS, "abs", 1e-12));
  }
  SECTION("G^BB_R fulfills the crossing relation in k-space") {
    opts.BB = RE;

    // Green's tensor with positive frequency and k-vector
    opts.omega = omega;
    opts.kvec(0) = phi;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    // Green's tensor with negative frequency and k-vector
    opts.omega = -omega;
    opts.kvec(0) = phi + M_PI;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        RHS(i, j) += (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_R(kappa_double, &opts);
        RHS(i, j) += I * (2. * M_PI) *
                     greens_tensor.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    REQUIRE(approx_equal(conj(LHS), RHS, "abs", 1e-12));
  }
  opts.BB = IGNORE;
  
}


TEST_CASE("Test the sum of several Green's tensors",
          "[GreensTensorPlateMagnetic]") {
  // In general there a lot of variations of Green's tensors that one could some
  // up. For simplicity only the two cases
  // which are relevant for quantum friction are tested, if to plan to use
  // other sums, consider testing them as well.
  // Here we considered also the volume element from the integration.
  std::complex<double> I(0.0, 1.0);
  // Calculate computes all possible entries of the electric Green's tensor. On
  // the other hand integrand_2d_k_magnetic already sets all entries linear in
  // k_y to zero, as the integral from 0 to 2 Pi would vanish there phi has to
  // be set, such that sin(phi) = 0
  double phi = GENERATE(0., M_PI);
  double k = GENERATE(take(3, random(0., 1e2)));
  // At this point only the evanescent part is implemented, therefore
  // omega_pl_quad has to be smaller then k_quad
  double omega = GENERATE(take(3, random(-.8, .8)));
  double k_x = k * cos(phi);
  double k_y = k * sin(phi);
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  double kappa_double;
  std::complex<double> kappa, volume_element;
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic Greens(v, za, 0.1, &refl, delta_cut, rel_err);
  struct Options_GreensTensorMagnetic opts;
  opts.class_pt = &Greens;

  // First, the calculate_tensor operation is used to generate the
  // Green's tensor with fancy_I
  cx_mat::fixed<3, 3> Green(fill::zeros);
  cx_mat::fixed<3, 3> Green_fancy_I(fill::zeros);
  cx_mat::fixed<3, 3> Green_fancy_R(fill::zeros);

  opts.kvec(0) = k_x;
  opts.kvec(1) = k_y;
  opts.omega = omega + k_x * v;
  k = sqrt(k_x * k_x + k_y * k_y);
  kappa = sqrt(std::complex<double>(k * k - opts.omega * opts.omega, 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));
  volume_element = kappa * k / (k - cos(phi) * v * opts.omega);

  Greens.calculate_tensor(Green, opts);
  if (opts.omega < 0) {
    volume_element = conj(volume_element);
  }
  Green *= volume_element;
  Green_fancy_I = (Green - trans(Green)) / (2. * I);
  Green_fancy_R = (Green + trans(Green)) / (2.);

  cx_mat::fixed<3, 3> Sigma = {
      {0, 0, 0}, {sin(phi), -cos(phi), 0}, {I * kappa / k, 0, -cos(phi)}};
  Sigma *= k * v / (opts.omega);

  cx_mat::fixed<3, 3> LHS(fill::zeros);
  cx_mat::fixed<3, 3> RHS(fill::zeros);

  if (kappa.real() == 0.) {
    kappa_double = kappa.imag();
  } else {
    kappa_double = kappa.real();

  } // Here we considered also the volume element from the integration.

  opts.kvec(0) = phi;
  opts.omega = omega;

  SECTION("Sum of G^EE_I and G^BE_R") {
    opts.fancy_complex = IM;
    opts.BE = RE;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) +=
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) +=
            (2 * M_PI) * I * Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
       }
      }
    RHS = Green_fancy_I + Sigma * Green_fancy_R;
    REQUIRE(approx_equal(LHS, RHS, "abs", 10e-12));
  }
  SECTION("Sum of fancy_I of all Green's tensors") {
    opts.fancy_complex = IM;
    opts.BE = IM;
    opts.EB = IM;
    opts.BB = IM;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) +=
            (2 * M_PI) * Greens.integrand_2d_k_magnetic_R(kappa_double, &opts);
        LHS(i, j) +=
            (2 * M_PI) * I * Greens.integrand_2d_k_magnetic_I(kappa_double, &opts);
      }
    }
    RHS = Green_fancy_I + Sigma * Green_fancy_I + Green_fancy_I * trans(Sigma) +
          Sigma * Green_fancy_I * trans(Sigma);
    REQUIRE(approx_equal(LHS, RHS, "abs", 10e-12));
  }
}


TEST_CASE("Integrand_1d_k_magnetic gives non-vanishing contriubtions for all"
          "Green's tensors ","GreensTensorPlateMagnetic")
{
  double omega = GENERATE(take(3,random(-1e2,1e2)));
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic Greens(v, za, 0.1, &refl, delta_cut, rel_err);
  struct Options_GreensTensorMagnetic opts;

  std::complex<double> I (0.,1.);

  auto phi = GENERATE(take(3,random(0.,2.*M_PI)));
  opts.class_pt = &Greens;
  opts.omega = omega;

  cx_mat::fixed<3, 3> LHS(fill::zeros);
  cx_mat::fixed<3, 3> RHS(fill::zeros);

  SECTION("G^EE_I") {
    opts.fancy_complex = IM;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) = Greens.integrand_1d_k_magnetic_R(phi, &opts);
        LHS(i, j) += I * Greens.integrand_1d_k_magnetic_I(phi, &opts);
      }
    }
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
  SECTION("G^EE_R") {
    opts.fancy_complex = RE;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) = Greens.integrand_1d_k_magnetic_R(phi, &opts);
        LHS(i, j) += I * Greens.integrand_1d_k_magnetic_I(phi, &opts);
      }
    }
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
  opts.fancy_complex = IGNORE;
  SECTION("G^BE_I") {
    opts.BE = IM;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) = Greens.integrand_1d_k_magnetic_R(phi, &opts);
        LHS(i, j) += I * Greens.integrand_1d_k_magnetic_I(phi, &opts);
      }
    }
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
  SECTION("G^BE_R") {
    opts.BE = RE;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) = Greens.integrand_1d_k_magnetic_R(phi, &opts);
        LHS(i, j) += I * Greens.integrand_1d_k_magnetic_I(phi, &opts);
      }
    }
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
  opts.BE = IGNORE;
  SECTION("G^EB_I") {
    opts.EB = IM;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) = Greens.integrand_1d_k_magnetic_R(phi, &opts);
        LHS(i, j) += I * Greens.integrand_1d_k_magnetic_I(phi, &opts);
      }
    }
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
  SECTION("G^EB_R") {
    opts.EB = RE;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) = Greens.integrand_1d_k_magnetic_R(phi, &opts);
        LHS(i, j) += I * Greens.integrand_1d_k_magnetic_I(phi, &opts);
      }
    }
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
  opts.BE = IGNORE;
  SECTION("G^BB_I") {
    opts.BB = IM;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) = Greens.integrand_1d_k_magnetic_R(phi, &opts);
        LHS(i, j) += I * Greens.integrand_1d_k_magnetic_I(phi, &opts);
      }
    }
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
  SECTION("G^BB_R") {
    opts.BB = RE;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;
        LHS(i, j) = Greens.integrand_1d_k_magnetic_R(phi, &opts);
        LHS(i, j) += I * Greens.integrand_1d_k_magnetic_I(phi, &opts);
      }
    }
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
}

TEST_CASE("integrand_2d_k for G^EE corresponds to GreensTensorPlate","[GreensTensorPlateMagnetic]")
{
  auto omega = GENERATE(231.21,0.9,-2.51, 9.73408);
  auto phi = GENERATE(0.,2.10922,1.3212);
  auto kappa_double = GENERATE(0.21, 21.1221, 201.12);

  std::complex<double> I (0.,1.);
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlate Greens(v, za, 0.1, &refl, delta_cut, rel_err);
  GreensTensorPlateMagnetic GreensM(v, za, 0.1, &refl, delta_cut, rel_err);
  struct Options_GreensTensorMagnetic opts;
  opts.omega = omega;
  opts.kvec(0) = phi;
  
  cx_mat::fixed<3,3> LHS(fill::zeros);
  cx_mat::fixed<3,3> RHS(fill::zeros);

  SECTION("IM") {
    opts.fancy_complex = IM;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;

	opts.class_pt = &GreensM;
	RHS(i,j) = GreensM.integrand_2d_k_magnetic_R(kappa_double,&opts);
	RHS(i,j) += I*GreensM.integrand_2d_k_magnetic_I(kappa_double, &opts);

	opts.class_pt = &Greens;
        LHS(i, j) = Greens.integrand_2d_k(kappa_double, &opts);
	if(i != j) LHS(i,j) *= I;
      }
    }

    REQUIRE(approx_equal(LHS, RHS, "abs", 1e-12));
  }
  
  SECTION("RE") {
    opts.fancy_complex = RE;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;

	opts.class_pt = &GreensM;
	RHS(i,j) = GreensM.integrand_2d_k_magnetic_R(kappa_double,&opts);
	RHS(i,j) += I*GreensM.integrand_2d_k_magnetic_I(kappa_double, &opts);

	opts.class_pt = &Greens;
        LHS(i, j) = Greens.integrand_2d_k(kappa_double, &opts);
	if(i != j) LHS(i,j) *= I;
      }
    }
    REQUIRE(approx_equal(LHS, RHS, "abs", 1e-11));
  }
  
}

TEST_CASE("integrand_1d_k of G^EE corresponds to GreensTensorPlate","[GreensTensorPlateMagnetic]")
{
  auto omega = GENERATE(take(1,random(-1e1,1e1)));
  auto phi = GENERATE(take(1,random(0.,2.*M_PI)));

  std::complex<double> I (0.,1.);
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlate Greens(v, za, 0.1, &refl, delta_cut, rel_err);
  GreensTensorPlateMagnetic GreensM(v, za, 0.1, &refl, delta_cut, rel_err);
  struct Options_GreensTensorMagnetic opts;
  opts.omega = omega;
  
  cx_mat::fixed<3,3> LHS(fill::zeros);
  cx_mat::fixed<3,3> RHS(fill::zeros);

  SECTION("IM") {
    opts.fancy_complex = IM;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;

	opts.class_pt = &GreensM;
	RHS(i,j) = GreensM.integrand_1d_k_magnetic_R(phi,&opts);
	RHS(i,j) += I*GreensM.integrand_1d_k_magnetic_I(phi, &opts);

	opts.class_pt = &Greens;
        LHS(i, j) = Greens.integrand_1d_k(phi, &opts);
	if(i != j) LHS(i,j) *= I;
      }
    }
    REQUIRE(approx_equal(LHS, RHS, "abs", 1e-12));
  }
  
  SECTION("RE") {
    opts.fancy_complex = RE;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        opts.indices(0) = i;
        opts.indices(1) = j;

	opts.class_pt = &GreensM;
	RHS(i,j) = GreensM.integrand_1d_k_magnetic_R(phi,&opts);
	RHS(i,j) += I*GreensM.integrand_1d_k_magnetic_I(phi, &opts);

	opts.class_pt = &Greens;
        LHS(i, j) = Greens.integrand_1d_k(phi, &opts);
	if(i != j) LHS(i,j) *= I;
      }
    }
    if(!approx_equal(LHS, RHS, "abs", 1e-9))
    {
      std::cout << "omega = " << omega << std::endl;
      std::cout << "phi = " << phi << std::endl;
      std::cout << LHS << RHS << std::endl;
      std::cout << LHS - RHS << std::endl;
    }
    REQUIRE(approx_equal(LHS, RHS, "abs", 1e-9));
  }
  
}
  


TEST_CASE("integrate_k gives non-vanishing contriubtions for all"
          "Green's tensors ","[GreensTensorPlateMagnetic]")
{
  double omega = GENERATE(take(3,random(-1e2,1e2)));
  double omega_p = 9;
  double gamma = .1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic Greens(v, za, 0.1, &refl, delta_cut, rel_err);
  struct Options_GreensTensorMagnetic opts;

  std::complex<double> I (0.,1.);

  opts.class_pt = &Greens;
  opts.omega = omega;

  cx_mat::fixed<3, 3> LHS(fill::zeros);
  cx_mat::fixed<3, 3> RHS(fill::zeros);

  SECTION("G^EE_I") {
    opts.fancy_complex = IM;
    Greens.integrate_k_magnetic(LHS, opts);
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
  SECTION("G^EE_R") {
    opts.fancy_complex = RE;
    Greens.integrate_k_magnetic(LHS, opts);
    REQUIRE(!approx_equal(LHS, RHS, "abs", 1e-12));
  }
}



TEST_CASE("Integrals of all the Green's tensor work properly",
          "[GreensTensorPlateMagnetic]") {
    double omega = GENERATE(take(3,random(0.,1e2)));
    double omega_p = 9;
    double gamma = .1;
    double v = 1e-2;
    double za = 0.1;
    double delta_cut = 100;
    vec::fixed<2> rel_err = {1E-8, 1E-6};
    PermittivityDrude perm(omega_p, gamma);
    ReflectionCoefficientsLocBulk refl(&perm);
    GreensTensorPlateMagnetic Greens(v, za, 0.1, &refl, delta_cut, rel_err);
    struct Options_GreensTensorMagnetic opts;
    opts.class_pt = &Greens;
    opts.omega = omega;

    cx_mat::fixed<3, 3> LHS(fill::zeros);
    cx_mat::fixed<3, 3> RHS(fill::zeros);
    SECTION("Integral over G^EE_I obeys the crossing relation") {
        opts.fancy_complex = IM;
        Greens.integrate_k_magnetic(LHS, opts);
        opts.omega = -omega;
        Greens.integrate_k_magnetic(RHS, opts);

	REQUIRE(approx_equal(LHS, -strans(RHS), "reldiff",10E-4));
    }

    
    SECTION("Integral over G^EE_R obeys the crossing relation") {
        opts.omega = omega;
        opts.fancy_complex = RE;
        Greens.integrate_k_magnetic(LHS, opts);

        opts.omega = -omega;
        Greens.integrate_k_magnetic(RHS, opts);

        REQUIRE(approx_equal(LHS, strans(RHS), "reldiff", 10E-4));
    }
    
    opts.fancy_complex = IGNORE;
    /*
    SECTION("Integral over G^BE_I obeys the crossing relation") {
        opts.omega = omega;
        opts.BE = IM;
        Greens.integrate_k_magnetic(LHS, opts);

        opts.omega = -omega;
        Greens.integrate_k_magnetic(RHS, opts);
	std::cout << LHS << -RHS << std::endl;
	std::cout << "Difference" << std::endl;
	std::cout << LHS + RHS << std::endl;

        REQUIRE(approx_equal(LHS, -RHS, "reldiff", 10E-4));
    }
 
    SECTION("Integral over G^BE_R obeys the crossing relation") {
        opts.omega = omega;
        opts.BE = RE;
        Greens.integrate_k_magnetic(LHS, opts);

        opts.omega = -omega;
        Greens.integrate_k_magnetic(RHS, opts);

        REQUIRE(approx_equal(LHS, RHS, "reldiff", 10E-4));
    }
    opts.BE = IGNORE;

  SECTION("Integral over G^EB_I obeys the crossing relation") {
    opts.omega = omega;
    opts.EB = IM;
    Greens.integrate_k_magnetic(LHS, opts);

    opts.omega = -omega;
    Greens.integrate_k_magnetic(RHS, opts);

    REQUIRE(approx_equal(LHS, -RHS, "reldiff", 10E-4));
  }

  SECTION("Integral over G^EB_R obeys the crossing relation") {
    opts.omega = omega;
    opts.EB = RE;
    Greens.integrate_k_magnetic(LHS, opts);

    opts.omega = -omega;
    Greens.integrate_k_magnetic(RHS, opts);

    REQUIRE(approx_equal(LHS, RHS, "reldiff", 10E-4));
  }
  SECTION("Integral over G^BB_I obeys the crossing relation") {
    opts.omega = omega;
    opts.BB = IM;
    Greens.integrate_k_magnetic(LHS, opts);

    opts.omega = -omega;
    Greens.integrate_k_magnetic(RHS, opts);

    REQUIRE(approx_equal(LHS, -RHS, "reldiff", 10E-4));
  }

  SECTION("Integral over G^BB_R obeys the crossing relation") {
    opts.omega = omega;
    opts.BB = RE;
    Greens.integrate_k_magnetic(LHS, opts);

    opts.omega = -omega;
    Greens.integrate_k_magnetic(RHS, opts);

    REQUIRE(approx_equal(LHS, RHS, "reldiff", 10E-4));
  }
  */

}


