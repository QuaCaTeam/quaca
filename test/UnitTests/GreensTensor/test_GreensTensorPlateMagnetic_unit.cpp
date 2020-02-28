#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>
#include <iomanip> // std::setprecision

TEST_CASE("integrand_2d_k_magnetic returns the correct Green's tensors",
          "[GreensTensorPlate]") {
  // Here we considered also the volume element from the integration.
  std::complex<double> I(0.0, 1.0);
  auto omega = GENERATE(take(5, random(-1e2, 1e2)));
  auto k = GENERATE(take(5, random(0., 1e2)));
  auto phi = GENERATE(take(5, random(0., 2*M_PI)));
  double omega_p = 9;
  double gamma = 0.1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  double kappa_double;
  double k_x, k_y;
  std::complex<double> kappa, volume_element;
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic greens_tensor(v, za, 0.1, &refl, delta_cut, rel_err);
  struct Options_GreensTensorMagnetic opts;
  opts.class_pt = &greens_tensor;

  // First, the calculate_tensor operation is used to generate the
  // Green's tensor with fancy_I
  cx_mat::fixed<3, 3> Green(fill::zeros);
  cx_mat::fixed<3, 3> Green_fancy_I(fill::zeros);
  cx_mat::fixed<3,3> Green_fancy_R(fill::zeros);

  opts.kvec(0) = phi;
  opts.kvec(1) = k;
  opts.omega = omega + k_x * v;
  k_x = cos(phi)*k;
  k_y = sin(phi)*k;
  kappa = sqrt(std::complex<double>(k * k - opts.omega * opts.omega, 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));
  //Volume element of the dk^2 integration
  volume_element = kappa * k / (k - cos(phi) * v * opts.omega);

  greens_tensor.calculate_tensor(Green, opts);
  if (opts.omega < 0) {
    volume_element = conj(volume_element);
  }
  Green *= volume_element;
  Green_fancy_I = (Green - trans(Green)) / (2. * I);
  Green_fancy_R = (Green + trans(Green)) / 2.;
  cx_mat::fixed<3,3> magnetic_tensor = {{0,0,0},
          {sin(phi),-cos(phi),0},
          {I*kappa/k,0,-cos(phi)}};
  magnetic_tensor *= k*v/omega;

  opts.kvec(0) = kappa_double;
  opts.kvec(1) = phi;

  SECTION("Test BE: RE and IM")
  {
      opts.BE = RE;
      //for loop to cycle through all elements
      for(int i = 0; i < 3; i++)
      {
          for(int j = 0; j < 3; j++)
          {
              opts.indices(0) = i;
              opts.indices(1) = j;
              greens_tensor.integrand_2d_k_magnetic(kappa_double, &opts);
          }
      }
      REQUIRE(approx_equal(Green, magnetic_tensor*Green_fancy_R, "reldiff", 1e-8));
  }

};

