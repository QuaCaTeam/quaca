#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Polarizability with bath can be constructed in different ways", "[PolarizabilityBath]")
{
  SECTION("Constructor from .ini file")
  {
    PolarizabilityBath pol("../data/test_files/PolarizabilityBath.ini");
    REQUIRE( pol.get_omega_a() == 1.3 );
    REQUIRE( pol.get_alpha_zero() == 6E-9 );

    // test if we read memory kernel correctly
    std::complex<double> test = pol.get_mu(3.0);
    REQUIRE( test.real() == 0.69420 );
  };

  SECTION("Constructor with direct input")
  {
    OhmicMemoryKernel mu(0.69420);
    PolarizabilityBath pol(1.3, 6E-9, &mu, NULL);
    REQUIRE( pol.get_omega_a() == 1.3 );
    REQUIRE( pol.get_alpha_zero() == 6E-9 );

    std::complex<double> test = pol.get_mu(3.0);
    REQUIRE( test.real() == 0.69420 );
  };
};

TEST_CASE("Check if integrand works", "[PolarizabilityBath]")
{
  // define greens tensor
  double v = 0.1;
  double beta = 10;
  GreensTensorVacuum greens(v, beta);

  // define polarizability
  double omega_a = 3.0;
  double alpha_zero = 2.4;
  double gamma = 2.0;
  OhmicMemoryKernel mu(gamma);
  PolarizabilityBath pol(omega_a, alpha_zero, &mu, &greens);

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

TEST_CASE("Test integration for omega_cut much smaller than omega_a", "[PolarizabilityNoBath]")
{
  // define greens tensor
  auto v = GENERATE(take(5, random(0.0, 1.0)));
  auto beta = GENERATE(take(5, random(0.0, 1e4)));
  GreensTensorVacuum greens(v, beta);

  // define polarizability
  auto omega_a = GENERATE(take(5, random(1.0, 1e3)));
  auto alpha_zero = GENERATE(take(5, random(0.0, 0.1)));
  auto gamma = GENERATE(take(5, random(0.0, 1e2)));
  OhmicMemoryKernel mu(gamma);
  PolarizabilityBath pol(omega_a, alpha_zero, &mu, &greens);

  Options_Polarizability opts;
  opts.fancy_I = true;
  opts.class_pt = &pol;

  double omega_min = 0.0;
  double fact = 1e-5;
  double omega_max = fact*omega_a;
  double relerr = 1e-12;
  double abserr = 0;

  double result, asymp; // double for result and asymptotic value

  // loop over indices
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      opts.indices(0) = i;
      opts.indices(1) = j;
      result = pol.integrate_omega(opts, omega_min, omega_max, relerr, abserr);

      // check diagonal entries
      if (i == j)
      {
        if (i == 0)
        {
          asymp = alpha_zero*alpha_zero*pow(omega_max,4)/2.0*1.0/(3*(1.0 - v*v)*(1.0 - v*v)) + alpha_zero*gamma*fact*fact/(2.0);
          REQUIRE(Approx(result).margin(relerr) == asymp);
        }
        else
        {
          asymp = alpha_zero*alpha_zero*pow(omega_max,4)/2.0*(1.0+v*v)/(3*pow((1.0 - v*v),3)) + alpha_zero*gamma*fact*fact/(2.0);
          REQUIRE(Approx(result).margin(relerr) == asymp);
        };
      }
      else
      {
        REQUIRE(result == 0); // off-diagonals are zero
      };

    };
  };

};


TEST_CASE("Test integration for omega_cut much larger than omega_a", "[PolarizabilityBath]")
{
  // define greens tensor
  double v = 1e-1;
  double beta = 1e3;
  GreensTensorVacuum greens(v, beta);

  // define polarizability
  double omega_a = 4.0;
  double alpha_zero = 1e-4;
  double gamma = 1.0;
  OhmicMemoryKernel mu(gamma);
  PolarizabilityBath pol(omega_a, alpha_zero, &mu, &greens);

  Options_Polarizability opts;
  opts.fancy_I = true;
  opts.indices(0) = 0;
  opts.indices(1) = 0;
  opts.class_pt = &pol;

  double omega_min = 0.0;
  double omega_max = omega_a*1e4;
  double relerr = 1e-15;
  double abserr = 0;

  double result = pol.integrate_omega(opts, omega_min, omega_a-1e-2, relerr, abserr);
  result += pol.integrate_omega(opts, omega_a-1e-2, omega_a+1e-2, relerr, abserr);
  result += pol.integrate_omega(opts, omega_a+1e-2, omega_max, relerr, abserr);
  double asymp = alpha_zero*omega_a*M_PI/2.0;

  //std::cout << result << std::endl;
  //std::cout << asymp << std::endl;

  REQUIRE(Approx(result).margin(relerr) == asymp);
};
