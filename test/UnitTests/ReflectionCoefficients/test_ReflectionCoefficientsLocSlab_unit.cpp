#include "Quaca.h"
#include "catch.hpp"
#include <complex>

TEST_CASE("ReflectionCoefficientsLocSlab constructors work as expected",
          "[ReflectionCoefficientsLocSlab]") {

  SECTION("Direct constructor") {
    double omega_p = 9;
    double omega = 1;
    double gamma = 3.5E-2;
    double thickness = 0.05;
    auto perm = std::make_shared<PermittivityDrude>(gamma, omega_p);
    ReflectionCoefficientsLocSlab RefC(perm, thickness);

    REQUIRE(RefC.get_epsilon(omega).real() == perm->calculate(omega).real());
    REQUIRE(RefC.get_thickness() == thickness);
  }

  SECTION("json file constructor") {
    double omega = 1;
    ReflectionCoefficientsLocSlab RefC(
        "../data/test_files/ReflectionLocalSlab.json");
    auto perm = std::make_shared<PermittivityDrude>(
        "../data/test_files/ReflectionLocalSlab.json");

    REQUIRE(RefC.get_epsilon(omega).real() == perm->calculate(omega).real());
    REQUIRE(RefC.get_thickness() == 0.05);
  }
}

TEST_CASE("For slab thickness*kappa >> 1, the bulk coefficient is retrieved",
          "[ReflectionCoefficientsLocSlab]") {
  auto omega = GENERATE(take(5, random(-1e-2, 1e-2)));

  double omega_p = 9;
  double gamma = 3.5E-2;
  auto perm = std::make_shared<PermittivityDrude>(gamma, omega_p);

  double thickness = 500;

  ReflectionCoefficientsLocSlab RefSlab(perm, thickness);
  ReflectionCoefficientsLocBulk RefBulk(perm);

  std::complex<double> kappa = 10.;

  std::complex<double> rp_slab, rs_slab;
  RefSlab.calculate(omega, kappa, rp_slab, rs_slab);
  std::complex<double> rp_bulk, rs_bulk;
  RefBulk.calculate(omega, kappa, rp_bulk, rs_bulk);

  REQUIRE(Approx(rp_slab.real()).epsilon(1E-4) == rp_bulk.real());
  REQUIRE(Approx(rp_slab.imag()).epsilon(1E-4) == rp_bulk.imag());
  REQUIRE(Approx(rs_slab.real()).epsilon(1E-4) == rs_bulk.real());
  REQUIRE(Approx(rs_slab.imag()).epsilon(1E-4) == rs_bulk.imag());
}

TEST_CASE("r_p and r_s of the slab configuration obey the crossing relation",
          "[ReflectionCoefficientsLocSlab]") {
  auto omega = GENERATE(-5.32,-0.12,1.1,8.6);

  auto kappa_double = GENERATE(-9.54,-1.32,3.2,8.9);
  std::complex<double> kappa;
  if (kappa_double < 0.) {
    kappa = std::complex<double>(0., -std::abs(kappa_double));
  } else {
    kappa = std::complex<double>(kappa_double, 0.);
  }

  ReflectionCoefficientsLocSlab RefC(
      "../data/test_files/ReflectionLocalSlab.json");

  std::complex<double> rp_lhs, rs_lhs;
  RefC.calculate(omega, kappa, rp_lhs, rs_lhs);

  std::complex<double> rp_rhs, rs_rhs;
  RefC.calculate(-omega, kappa, rp_rhs, rs_rhs);

  REQUIRE(rp_lhs.real() == rp_rhs.real());
  REQUIRE(rs_lhs.imag() == -rs_rhs.imag());
}
