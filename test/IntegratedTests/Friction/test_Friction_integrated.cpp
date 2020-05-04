#include "Quaca.h"
#include "catch.hpp"
#include <iostream>

TEST_CASE("Analytical results with vacuum Green's tensor gets reproduced",
          "[Friction]") {
  // Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  double omega_a = .3;
  double alpha_zero = 6e-9;

  double beta = 1e-1;
  double v = 1e-5;
  double analytical_result = -2. / 12. * v * alpha_zero * pow(omega_a, 6) *
                             beta / pow(sinh(omega_a * beta / 2.), 2);
  // double analytical_result = (3./2.)*alpha_zero*pow(omega_a,2)/beta;

  double relerr_omega = 1e-6;

  double relerr_k = 1E-9;
  auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);
  auto alpha = std::make_shared<Polarizability>(omega_a, alpha_zero, greens);
  auto powerspectrum = std::make_shared<PowerSpectrum>(greens, alpha);
  Friction quant_fric(greens, alpha, powerspectrum, relerr_omega);

  double num_result = quant_fric.calculate(NON_LTE_ONLY);
  REQUIRE(Approx(num_result).epsilon(1e-4) == analytical_result);
}

TEST_CASE("Analytical results with scattered Green's tensor gets reproduced",
          "[Friction]") {
  // Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  double omega_a = 1.3;
  double alpha_zero = 6e-9;
  double za = 0.01;
  double omega_p = 9.;
  double gamma = 0.1;
  double rho;
  rho = gamma * M_PI * 4. / pow(omega_p, 2);
  double beta = 1e6;
  double v = 1e-4;
  double delta_cut = 30;
  double analytical_result = -(63. - 45.) * pow(alpha_zero * rho, 2) *
                             pow(v / M_PI, 3) / pow(2 * za, 10);
  analytical_result += -(6. - 3.) * pow(alpha_zero * rho / beta, 2) *
                       (v / M_PI) / pow(2 * za, 8);

  vec::fixed<2> rel_err = {1E-6, 1E-4};
  double relerr_omega = 1e-2;

  auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
  auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);
  auto greens = std::make_shared<GreensTensorPlate>(v, beta, za, refl,
                                                    delta_cut, rel_err);
  auto alpha = std::make_shared<Polarizability>(omega_a, alpha_zero, greens);
  auto powerspectrum = std::make_shared<PowerSpectrum>(greens, alpha);
  Friction quant_fric(greens, alpha, powerspectrum, relerr_omega);

  double num_result = quant_fric.calculate(NON_LTE_ONLY);
  std::cout << "ana=" << analytical_result << std::endl;
  std::cout << "num=" << num_result << std::endl;
  std::cout << "num/ana=" << num_result / analytical_result << std::endl;
  REQUIRE(Approx(num_result).epsilon(1e-2) == analytical_result);
}
