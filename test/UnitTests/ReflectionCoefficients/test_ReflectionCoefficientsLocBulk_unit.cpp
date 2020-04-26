#include "Quaca.h"
#include "catch.hpp"
#include <complex>

TEST_CASE("ReflectionCoefficientsLocBulk constructors work as expected",
          "[ReflectionCoefficientsLocBulk]") {

  SECTION("Direct constructor") {
    double omega_p = 9;
    double omega = 1;
    double gamma = 3.5E-2;
    auto perm = std::make_shared<PermittivityDrude>(gamma, omega_p);
    ReflectionCoefficientsLocBulk RefC(perm);

    REQUIRE(Approx(RefC.get_epsilon(omega).real()).epsilon(1E-6) ==
            perm->calculate(omega).real());
  }

  SECTION("Constructor with json file works") {
    double omega = 1;
    ReflectionCoefficientsLocBulk RefC(
        "../data/test_files/GreensTensorPlate.json");
    auto perm = std::make_shared<PermittivityDrude>(
        "../data/test_files/GreensTensorPlate.json");

    REQUIRE(Approx(RefC.get_epsilon(omega).real()).epsilon(1E-6) ==
            perm->calculate(omega).real());
  }
}

TEST_CASE("ReflectionCoefficientsLocBulk reproduces evanescent limit",
          "[ReflectionCoefficientsLocBulk]") {
  auto omega = GENERATE(take(5, random(-1e-2, 1e-2)));
  std::complex<double> kappa = 100.;
  ReflectionCoefficientsLocBulk RefC(
      "../data/test_files/GreensTensorPlate.json");
  auto perm = std::make_shared<PermittivityDrude>(
      "../data/test_files/GreensTensorPlate.json");

  double k_quad = std::real(kappa * kappa) + omega * omega;

  std::complex<double> rp, rs;
  RefC.calculate(omega, kappa, rp, rs);

  std::complex<double> rp_app =
      (perm->calculate(omega) - 1.) / (perm->calculate(omega) + 1.);
  std::complex<double> rs_app =
      (perm->calculate(omega) - 1.) * omega * omega / (4 * k_quad);
  rs_app += (pow(perm->calculate(omega), 2) - 1.) * pow(omega, 4) /
            (8 * k_quad * k_quad);
  REQUIRE(Approx(rp.real()).epsilon(1E-4) == rp_app.real());
  REQUIRE(Approx(rp.imag()).epsilon(1E-4) == rp_app.imag());
  REQUIRE(Approx(rs.real()).epsilon(1E-4) == rs_app.real());
  REQUIRE(Approx(rs.imag()).epsilon(1E-4) == rs_app.imag());
}

TEST_CASE("ReflectionCoefficientsLocBulk reproduces propagation limit",
          "[ReflectionCoefficientsLocBulk]") {
  //    auto omega = GENERATE(take(5, random(1e1, 1e2)));
  double omega = 100;
  double k = 0.1;

  ReflectionCoefficientsLocBulk RefC(
      "../data/test_files/GreensTensorPlate.json");
  auto perm = std::make_shared<PermittivityDrude>(
      "../data/test_files/GreensTensorPlate.json");

  std::complex<double> eps = perm->calculate(omega);

  std::complex<double> rp, rs;
  std::complex<double> kappa =
      std::complex<double>(0., -sqrt(omega * omega - k * k));
  RefC.calculate(omega, kappa, rp, rs);

  std::complex<double> rs_app = (1. - sqrt(eps)) / (1. + sqrt(eps));
  std::complex<double> rp_app = (eps - sqrt(eps)) / (eps + sqrt(eps));

  REQUIRE(Approx(rp.real()).epsilon(1E-4) == rp_app.real());
  REQUIRE(Approx(rp.imag()).epsilon(1E-4) == rp_app.imag());
  REQUIRE(Approx(rs.real()).epsilon(1E-4) == rs_app.real());
  REQUIRE(Approx(rs.imag()).epsilon(1E-4) == rs_app.imag());
}

TEST_CASE("r_p and r_s of bulk obey the crossing relation",
          "[ReflectionCoefficientsLocBulk]") {
  auto omega = GENERATE(take(5, random(-1e2, 1e2)));
  auto kappa_double = GENERATE(take(5, random(-1e2, 1e2)));

  std::complex<double> kappa;
  if (kappa_double < 0.) {
    kappa = std::complex<double>(0., -std::abs(kappa_double));
  } else {
    kappa = std::complex<double>(kappa_double, 0.);
  }

  ReflectionCoefficientsLocBulk RefC(
      "../data/test_files/GreensTensorPlate.json");
  auto perm = std::make_shared<PermittivityDrude>(
      "../data/test_files/GreensTensorPlate.json");

  std::complex<double> rp_lhs, rs_lhs;
  RefC.calculate(omega, kappa, rp_lhs, rs_lhs);

  std::complex<double> rp_rhs, rs_rhs;
  RefC.calculate(-omega, kappa, rp_rhs, rs_rhs);

  REQUIRE(rp_lhs.real() == rp_rhs.real());
  REQUIRE(rs_lhs.imag() == -rs_rhs.imag());
}
