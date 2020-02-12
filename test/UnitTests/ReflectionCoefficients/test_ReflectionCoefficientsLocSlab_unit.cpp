#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>

TEST_CASE("ReflectionCoefficientsLocSlab constructors work as expected",
          "[ReflectionCoefficientsLocSlab]") {

  SECTION("Direct constructor") {
    double omega_p = 9;
    double omega = 1;
    double gamma = 3.5E-2;
    double thickness = 0.05;
    Permittivity *perm = new PermittivityDrude(gamma, omega_p);
    ReflectionCoefficientsLocSlab RefC(perm, thickness);

    REQUIRE(RefC.get_epsilon(omega).real() == perm->epsilon(omega).real());
    REQUIRE(RefC.get_thickness() == thickness);
  };

  SECTION("ini file constructor") {
    double omega = 1;
    ReflectionCoefficientsLocSlab RefC(
        "../data/test_files/GreensTensorSlab.ini");
    Permittivity *perm =
        new PermittivityDrude("../data/test_files/GreensTensorSlab.ini");

    REQUIRE(RefC.get_epsilon(omega).real() == perm->epsilon(omega).real());
    REQUIRE(RefC.get_thickness() == 0.05);
  };
};

TEST_CASE("For slab thickness*kappa >> 1, the bulk coefficient is retrieved",
          "[ReflectionCoefficientsLocSlab]") {
  auto omega = GENERATE(take(5, random(-1e-2, 1e-2)));
  std::complex<double> kappa = 10.;
  std::complex<double> rp_bulk, rs_bulk, rp_slab, rs_slab;
  double omega_p = 9;
  double gamma = 3.5E-2;
  double thickness = 500;
  Permittivity *perm = new PermittivityDrude(gamma, omega_p);
  ReflectionCoefficientsLocSlab RefSlab(perm, thickness);
  ReflectionCoefficientsLocBulk RefBulk(perm);
  double k_quad = std::real(kappa * kappa) + omega * omega;

  RefSlab.ref(rp_slab, rs_slab, omega, kappa);
  RefBulk.ref(rp_bulk, rs_bulk, omega, kappa);
  REQUIRE(Approx(rp_slab.real()).epsilon(1E-4) == rp_bulk.real());
  REQUIRE(Approx(rp_slab.imag()).epsilon(1E-4) == rp_bulk.imag());
  REQUIRE(Approx(rs_slab.real()).epsilon(1E-4) == rs_bulk.real());
  REQUIRE(Approx(rs_slab.imag()).epsilon(1E-4) == rs_bulk.imag());
};

TEST_CASE("r_p and r_s of the slab configuration obey the crossing relation",
          "[ReflectionCoefficientsLocSlab]") {
  auto omega = GENERATE(take(5, random(-1e2, 1e2)));
  auto kappa_double = GENERATE(take(5, random(-1e2, 1e2)));
  std::complex<double> kappa;
  if (kappa_double < 0.) {
    kappa = std::complex<double>(0., -std::abs(kappa_double));
  } else {
    kappa = std::complex<double>(kappa_double, 0.);
  };
  std::complex<double> rp_lhs, rs_lhs, rp_rhs, rs_rhs;
  ReflectionCoefficientsLocSlab RefC("../data/test_files/GreensTensorSlab.ini");
  Permittivity *perm =
      new PermittivityDrude("../data/test_files/GreensTensorSlab.ini");

  RefC.ref(rp_lhs, rs_lhs, omega, kappa);
  RefC.ref(rp_rhs, rs_rhs, -omega, kappa);

  REQUIRE(rp_lhs.real() == rp_rhs.real());
  REQUIRE(rs_lhs.imag() == -rs_rhs.imag());
};
