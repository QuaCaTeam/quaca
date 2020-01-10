#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Check if integrand works for no bath", "[PolarizabilityNoBath]")
{
  // define greens tensor
  double v = 0.1;
  double beta = 10;
  GreensTensorVacuum greens(v, beta);

  // define polarizability
  double omega_a = 3.0;
  double alpha_zero = 2.4;
  PolarizabilityNoBath pol(omega_a, alpha_zero, &greens);

  // frequency to evaluate
  double omega = 3.0;

  // define options struct for integrand
  Options_Polarizability opts;
  opts.class_pt = &pol;
  opts.fancy_I = true;
  opts.indices(0) = 0;
  opts.indices(1) = 0;

  // calculate as a reference the normal way
  cx_mat::fixed<3,3> alpha;
  opts.omega = omega;
  pol.calculate_tensor(alpha, opts);


  // calculation with calculate_tensor give the same result as integrand

  // loop over all indices
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      opts.indices(0) = i;
      opts.indices(1) = j;
      pol.calculate_tensor(alpha, opts);
      REQUIRE(pol.integrand_omega(omega, &opts) == alpha(opts.indices(0), opts.indices(1)));
    }
  }

};

TEST_CASE("Test integration for omega_cut much smaller than omega_a for no bath", "[PolarizabilityNoBath]")
{
  // define greens tensor
  double v = 0.4;
  double beta = 1e2;
  GreensTensorVacuum greens(v, beta);

  // define polarizability
  double omega_a = 4.0;
  double alpha_zero = 2.4;
  PolarizabilityNoBath pol(omega_a, alpha_zero, &greens);

  Options_Polarizability opts;
  opts.fancy_I = true;
  opts.indices(0) = 0;
  opts.indices(1) = 0;
  opts.class_pt = &pol;

  double omega_min = 0.0;
  double omega_max = 1e-2;
  double relerr = 1e-12;
  double abserr = 0;

  double result = pol.integrate_omega(opts, omega_min, omega_max, relerr, abserr);
  double asymp = alpha_zero*alpha_zero*pow(omega_max,4)/2.0*1.0/(3*(1.0 - v*v)*(1.0 - v*v));

  //std::cout << result << std::endl;
  //std::cout << asymp << std::endl;

  REQUIRE(Approx(result).margin(relerr) == asymp);
};


TEST_CASE("Test integration for omega_cut much larger than omega_a for no bath", "[PolarizabilityNoBath]")
{
  // define greens tensor
  double v = 0.4;
  double beta = 1e-1;
  GreensTensorVacuum greens(v, beta);

  // define polarizability
  double omega_a = 4.0;
  double alpha_zero = 1e-6;
  PolarizabilityNoBath pol(omega_a, alpha_zero, &greens);

  Options_Polarizability opts;
  opts.fancy_I = true;
  opts.indices(0) = 0;
  opts.indices(1) = 0;
  opts.class_pt = &pol;

  double omega_min = 0.0;
  double omega_max = omega_a*1e4;
  double relerr = 1e-15;
  double abserr = 0;

  double result = pol.integrate_omega(opts, omega_min, omega_a, relerr, abserr);
  result += pol.integrate_omega(opts, omega_a, omega_max, relerr, abserr);
  double asymp = alpha_zero*omega_a*M_PI/2.0;

  //std::cout << result << std::endl;
  //std::cout << asymp << std::endl;

  REQUIRE(Approx(result).margin(1e-7) == asymp);
};
