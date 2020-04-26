#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>

TEST_CASE("Plate Green's tensor constructors work as expected",
          "[GreensTensorPlate]") {

  SECTION("Direct constructor") {

    double omega_p = 9;
    double gamma = 0.1;
    auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
    auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);

    double v = 1E-5;
    double za = 0.1;
    double beta = 1E4;
    double delta_cut = 20;
    vec::fixed<2> rel_err = {1E-8, 1E-6};
    GreensTensorPlate Greens(v, beta, za, refl, delta_cut, rel_err);

    REQUIRE(Greens.get_za() == za);
    REQUIRE(Greens.get_delta_cut() == delta_cut);
    REQUIRE(Greens.get_rel_err_0() == rel_err(0));
    REQUIRE(Greens.get_rel_err_1() == rel_err(1));

    double omega = 1;
    double k = 10;
    std::complex<double> kappa;
    if (k < omega) {
      kappa = std::complex<double>(0., -sqrt(omega * omega - k * k));
    } else {
      kappa = std::complex<double>(sqrt(k * k - omega * omega), 0.);
    }

    std::complex<double> rp, rs;
    refl->calculate(omega, kappa, rp, rs);

    REQUIRE(real(Greens.get_r_s(omega, k)) == rs.real());
    REQUIRE(real(Greens.get_r_p(omega, k)) == rp.real());
    REQUIRE(imag(Greens.get_r_s(omega, k)) == rs.imag());
    REQUIRE(imag(Greens.get_r_p(omega, k)) == rp.imag());
  }

  SECTION("json file constructor") {
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.json");

    REQUIRE(Greens.get_za() == 0.1);
    REQUIRE(Greens.get_delta_cut() == 20);
    REQUIRE(Greens.get_rel_err_0() == 1E-8);
    REQUIRE(Greens.get_rel_err_1() == 1E-6);

    std::complex<double> rp, rs;
    double omega = 1;
    double k = 10;
    std::complex<double> kappa;
    if (k < omega) {
      kappa = std::complex<double>(0., -sqrt(omega * omega - k * k));
    } else {
      kappa = std::complex<double>(sqrt(k * k - omega * omega), 0.);
    }

    ReflectionCoefficientsLocBulk refl(
        "../data/test_files/GreensTensorPlate.json");
    refl.calculate(omega, kappa, rp, rs);

    REQUIRE(real(Greens.get_r_s(omega, k)) == rs.real());
    REQUIRE(real(Greens.get_r_p(omega, k)) == rp.real());
    REQUIRE(imag(Greens.get_r_s(omega, k)) == rs.imag());
    REQUIRE(imag(Greens.get_r_p(omega, k)) == rp.imag());
  }
}

TEST_CASE("The operations calculate_tensor and the integrand_2d_k coincide",
          "[GreensTensorPlate]") {
  // Here we considered also the volume element from the integration.
  std::complex<double> I(0.0, 1.0);
  auto omega = GENERATE(take(5, random(-1e2, 1e2)));
  auto k_x = GENERATE(take(5, random(0., 1e2)));
  auto k_y = GENERATE(take(5, random(0., 1e2)));
  double omega_p = 9;
  double gamma = 0.1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};

  auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
  auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);
  GreensTensorPlate Greens(v, za, 0.1, refl, delta_cut, rel_err);

  // First, the calculate_tensor operation is used to generate the
  // Green's tensor with fancy_I
  cx_mat::fixed<3, 3> Green(fill::zeros);
  cx_mat::fixed<3, 3> Green_fancy_I_ct(fill::zeros);

  vec::fixed<2> kvec;
  kvec(0) = k_x;
  kvec(1) = k_y;
  omega = omega + k_x * v;
  double k = sqrt(k_x * k_x + k_y * k_y);
  double cos_phi = k_x / k;
  std::complex<double> kappa =
      sqrt(std::complex<double>(k * k - omega * omega, 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));
  double volume_element = std::abs(kappa) * k / (k - cos_phi * v * omega);

  Greens.calculate_tensor(omega, kvec, Green);

  Green *= volume_element;
  Green_fancy_I_ct = (Green - trans(Green)) / (2. * I);

  // Second, the integrand_2d_k operation is used
  cx_mat::fixed<3, 3> Green_fancy_I_ik2d(fill::zeros);
  double kappa_double;
  if (kappa.real() == 0.) {
    kappa_double = kappa.imag();
  } else {
    kappa_double = kappa.real();
  }

  kvec(0) = acos(k_x / k);
  omega = omega - k_x * v;

  // loop over all indices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Green_fancy_I_ik2d(i, j) =
          (2 * M_PI) * Greens.integrand_2d_k(kappa_double, omega, kvec(0),
                                             {(double)i, (double)j}, IM, UNIT);
      if (i != j) {
        // As the prefactor I can not be evaluated in a purely real integration
        // routine, it was dropped in integrate_2d_k and has to be inserted here
        Green_fancy_I_ik2d(i, j) *= I;
      }
    }
  }
  REQUIRE(approx_equal(Green_fancy_I_ct, Green_fancy_I_ik2d, "reldiff", 10E-4));
}

TEST_CASE("Plate Green's tensor fulfills physical relations",
          "[GreensTensorPlate]") {

  SECTION("Green's tensor obeys reciprocity") {
    auto omega = GENERATE(take(5, random(0., 1e2)));
    auto k_x = GENERATE(take(5, random(0., 1e2)));
    auto k_y = GENERATE(take(5, random(0., 1e2)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.json");

    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

    Greens.calculate_tensor(omega, {k_x, k_y}, Greens_lhs);
    Greens.calculate_tensor(omega, {-k_x, -k_y}, Greens_rhs);

    REQUIRE(approx_equal(Greens_lhs, strans(Greens_rhs), "reldiff", 10E-15));
  }

  SECTION("Green's tensor obeys reality condition") {
    auto omega = GENERATE(take(5, random(0., 1e2)));
    auto k_x = GENERATE(take(5, random(0., 1e2)));
    auto k_y = GENERATE(take(5, random(0., 1e2)));

    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.json");

    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

    Greens.calculate_tensor(omega, {k_x, k_y}, Greens_lhs);
    Greens.calculate_tensor(-omega, {k_x, k_y}, Greens_rhs);

    REQUIRE(approx_equal(Greens_lhs, trans(Greens_rhs), "reldiff", 10E-15));
  }
}

TEST_CASE("Integrated Green's tensor works properly", "[GreensTensorPlate]") {

  SECTION("Integral over Green_fancy_I obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e2)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.json");

    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

    Greens.integrate_k(omega, Greens_lhs, IM, UNIT);
    Greens.integrate_k(-omega, Greens_rhs, IM, UNIT);

    REQUIRE(approx_equal(Greens_lhs, -strans(Greens_rhs), "reldiff", 10E-4));
  }

  SECTION("Integral over Green_fancy_R obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e2)));

    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.json");

    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

    Greens.integrate_k(omega, Greens_lhs, RE, UNIT);
    Greens.integrate_k(-omega, Greens_rhs, RE, UNIT);

    REQUIRE(approx_equal(Greens_lhs, strans(Greens_rhs), "reldiff", 10E-4));
  }

  SECTION("Integral over Green_fancy_I_kv obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e2)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.json");

    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);

    Greens.integrate_k(omega, Greens_lhs, IM, KV);
    Greens.integrate_k(-omega, Greens_rhs, IM, KV);

    REQUIRE(approx_equal(Greens_lhs, strans(Greens_rhs), "reldiff", 10E-4));
  }

  SECTION("Integral over Green fancy_I, temp obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e-1)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.json");

    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs1(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs2(fill::zeros);

    Greens.integrate_k(-omega, Greens_lhs, IM, TEMP);
    Greens.integrate_k(omega, Greens_rhs1, IM, UNIT);
    Greens.integrate_k(omega, Greens_rhs2, IM, TEMP);
    Greens_rhs = -Greens_rhs1 + Greens_rhs2;
    REQUIRE(approx_equal(Greens_lhs, strans(Greens_rhs), "reldiff", 10E-4));
  }

  SECTION("Integral over Green_fancy_I_kv_temp obeys the crossing relation") {
    auto omega = GENERATE(take(1, random(0., 1e-1)));
    GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.json");

    cx_mat::fixed<3, 3> Greens_lhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs1(fill::zeros);
    cx_mat::fixed<3, 3> Greens_rhs2(fill::zeros);

    Greens.integrate_k(-omega, Greens_lhs, IM, KV_TEMP);
    Greens.integrate_k(omega, Greens_rhs1, IM, KV);
    Greens.integrate_k(omega, Greens_rhs2, IM, KV_TEMP);

    Greens_rhs = -Greens_rhs1 + Greens_rhs2;
    REQUIRE(approx_equal(Greens_lhs, -strans(Greens_rhs), "reldiff", 10E-4));
  }
}

TEST_CASE("Integrated Green's tensor matches asymptotes",
          "[GreensTensorPlate]") {
  SECTION("Low-frequency asymptote of fancy_I") {
    // \omega << \omega_p and \omega << v/z_a
    std::complex<double> I(0.0, 1.0);
    double omega_p = 9;
    double gamma = 0.1;
    double v = 1e-5;
    double za = 0.1;
    auto omega = GENERATE(take(10, random(-0.1 * 1e-6, 0.1 * 1e-6)));
    double delta_cut = 30;
    vec::fixed<2> rel_err = {1E-8, 1E-6};

    auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
    auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);
    GreensTensorPlate Greens(v, 0.1, za, refl, delta_cut, rel_err);

    cx_mat::fixed<3, 3> GT_Ana(fill::zeros);
    GT_Ana(0, 0) = 2 * omega * gamma / pow(omega_p, 2) / pow(2 * za, 3);
    GT_Ana(1, 1) = GT_Ana(0, 0);
    GT_Ana(2, 2) = 2. * GT_Ana(0, 0);
    GT_Ana(0, 2) = (2 * 3 * v * gamma / (pow(omega_p, 2) * pow(2 * za, 4))) / I;
    GT_Ana(2, 0) = -GT_Ana(0, 2);

    cx_mat::fixed<3, 3> GT_Num(fill::zeros);
    Greens.integrate_k(omega, GT_Num, IM, UNIT);

    REQUIRE(approx_equal(GT_Ana, GT_Num, "reldiff", 10E-4));
  }

  SECTION("Low-frequency and low temperature asymptote of fancy_I_temp") {
    // \omega << \omega_p and \omega << v/z_a
    std::complex<double> I(0.0, 1.0);
    double omega_p = 9;
    double gamma = 0.1;
    double v = 1e-5;
    double za = 0.1;
    double beta = 1e12;
    auto omega = GENERATE(take(1, random(0., 1e-6)));
    double delta_cut = 30;
    vec::fixed<2> rel_err = {1E-8, 1E-6};

    auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
    auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);
    GreensTensorPlate Greens(v, beta, za, refl, delta_cut, rel_err);

    cx_mat::fixed<3, 3> GT_Ana(fill::zeros);
    double eta = omega * 2 * za / v;
    double rho = gamma / pow(omega_p, 2);
    GT_Ana(0, 0) =
        (v * rho * 2. / (pow(2 * za, 4) * M_PI)) * (0.5 * M_PI * eta + 4.);
    GT_Ana(1, 1) =
        (v * rho * 2. / (pow(2 * za, 4) * M_PI)) * (0.5 * M_PI * eta + 2.);
    GT_Ana(2, 2) = (v * rho * 2. / (pow(2 * za, 4) * M_PI)) * (M_PI * eta + 6.);
    GT_Ana(2, 0) = I * (v * rho * 2. / (pow(2 * za, 4) * M_PI)) *
                   (3. / 2. * M_PI + 2. * eta);
    GT_Ana(0, 2) = -GT_Ana(2, 0);

    cx_mat::fixed<3, 3> GT_Num(fill::zeros);
    Greens.integrate_k(omega, GT_Num, IM, TEMP);

    REQUIRE(approx_equal(GT_Ana, GT_Num, "reldiff", 10E-4));
  }

  SECTION("Low-frequency and low temperature asymptote of fancy_I_temp_kv") {
    // \omega << \omega_p and \omega << v/z_a
    std::complex<double> I(0.0, 1.0);
    double omega_p = 9;
    double gamma = 0.1;
    double v = 1e-5;
    double za = 0.1;
    double beta = 1e12;
    auto omega = GENERATE(take(1, random(0., 1e-6)));
    double delta_cut = 30;
    vec::fixed<2> rel_err = {1E-8, 1E-6};

    auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
    auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);
    GreensTensorPlate Greens(v, beta, za, refl, delta_cut, rel_err);

    cx_mat::fixed<3, 3> GT_Ana(fill::zeros);
    double eta = omega * 2 * za / v;
    double rho = gamma / pow(omega_p, 2);
    GT_Ana(0, 0) =
        (v * rho * 2. / (pow(2 * za, 5) * M_PI)) * (0.5 * M_PI * 9 + 4. * eta);
    GT_Ana(1, 1) =
        (v * rho * 2. / (pow(2 * za, 5) * M_PI)) * (0.5 * M_PI * 3 + 2. * eta);
    GT_Ana(2, 2) = GT_Ana(0, 0) + GT_Ana(1, 1);
    GT_Ana(2, 0) = I * (v * rho * 2. / (pow(2 * za, 5) * M_PI)) *
                   (3. / 2. * M_PI * eta + 16.);
    GT_Ana(0, 2) = -GT_Ana(2, 0);

    cx_mat::fixed<3, 3> GT_Num(fill::zeros);
    Greens.integrate_k(omega, GT_Num, IM, KV_TEMP);

    REQUIRE(approx_equal(GT_Ana, GT_Num, "reldiff", 10E-4));
  }

  SECTION("Low-frequency and high temperature asymptote of fancy_I_temp") {
    // \omega << \omega_p and \omega << v/z_a
    std::complex<double> I(0.0, 1.0);
    double omega_p = 9;
    double gamma = 0.1;
    double v = 1e-5;
    double za = 0.1;
    double beta = 1e-1;
    auto omega = GENERATE(take(1, random(1e-8, 1e-7)));
    double delta_cut = 30;
    vec::fixed<2> rel_err = {1E-8, 1E-6};

    auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
    auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);
    GreensTensorPlate Greens(v, beta, za, refl, delta_cut, rel_err);

    cx_mat::fixed<3, 3> GT_lhs(fill::zeros);
    Greens.integrate_k(omega, GT_lhs, IM, TEMP);

    cx_mat::fixed<3, 3> GT_rhs(fill::zeros);
    GT_rhs(0, 0) = 2 * gamma / (pow(omega_p, 2) * pow(2 * za, 3) * beta);
    GT_rhs(1, 1) = GT_rhs(0, 0);
    GT_rhs(2, 2) = GT_rhs(1, 1) + GT_rhs(0, 0);
    GT_rhs(0, 2) =
        0.5 * (2 * 3 * v * gamma / (pow(omega_p, 2) * pow(2 * za, 4))) / I;
    GT_rhs(2, 0) = -GT_rhs(0, 2);

    REQUIRE(approx_equal(GT_lhs, GT_rhs, "reldiff", 10E-4));
  }
}
