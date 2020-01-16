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
  auto v = GENERATE(take(3, random(0.0, 1.0)));
  auto beta = GENERATE(take(3, random(0.0, 1e4)));
  GreensTensorVacuum greens(v, beta);

  // define polarizability
  auto omega_a = GENERATE(take(3, random(1e-1, 1e1)));
  auto alpha_zero = GENERATE(take(3, random(0.0, 0.1)));
  PolarizabilityNoBath pol(omega_a, alpha_zero, &greens);

  Options_Polarizability opts;
  opts.fancy_I = true;
  opts.class_pt = &pol;

  double omega_min = 0.0;
  double omega_max = 1e-3*omega_a; // omega_max much smaller than omega_a
  double relerr = 1e-13;
  double abserr = 0;

  double result, asymp; // double for result and asymptotic value
  double toterr;

  // loop over indices
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      opts.indices(0) = i;
      opts.indices(1) = j;
      result = pol.integrate_omega(opts, omega_min, omega_max, relerr, abserr);

      /*
      * error tolerance includes absolute error from integration, which is result*relerr
      * and error from series expansion in omega_a.
      * we roughly estimate the series error with (omega_max)^2
      */
      toterr = result*relerr + omega_max*omega_max;


      // check diagonal entries
      if (i == j)
      {
        if (i == 0)
        {
          asymp = alpha_zero*alpha_zero*pow(omega_max,4)/2.0*1.0/(3*(1.0 - v*v)*(1.0 - v*v));
          REQUIRE(Approx(result).margin(toterr) == asymp);
        }
        else
        {
          asymp = alpha_zero*alpha_zero*pow(omega_max,4)/2.0*(1.0+v*v)/(3*pow((1.0 - v*v),3));
          REQUIRE(Approx(result).margin(toterr) == asymp);
        };
      }
      else
      {
        REQUIRE(result == 0); // off-diagonals are zero
      };

    };
  };

};


TEST_CASE("Test integration for omega_cut much larger than omega_a for no bath", "[PolarizabilityNoBath]")
{
  // define greens tensor
  auto v = GENERATE(take(3, random(0.0, 1.0)));
  double beta = 1e5;
  GreensTensorVacuum greens(v, beta);

  // define polarizability
  auto omega_a = GENERATE(take(3, random(1e-1, 1e1)));
  double alpha_zero = 1e-10;
  PolarizabilityNoBath pol(omega_a, alpha_zero, &greens);

  Options_Polarizability opts;
  opts.fancy_I = true;
  opts.class_pt = &pol;

  double omega_min = 0.0;
  double omega_max = omega_a*1e5;
  double relerr = 1e-13;
  double abserr = 0.;

  double result;
  double asymp = alpha_zero*omega_a*M_PI/2.0;

  double toterr;


  // loop over indices
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      opts.indices(0) = i;
      opts.indices(1) = j;
      result = pol.integrate_omega(opts, omega_min, omega_a-1e-3, relerr, abserr);
      result += pol.integrate_omega(opts, omega_a-1e-3, omega_a+1e-3, relerr, abserr);
      result += pol.integrate_omega(opts, omega_a+1e-3, omega_max, relerr, abserr);

      /*
      * error tolerance includes absolute error from integration, which is result*relerr
      */
      toterr = result*relerr + alpha_zero*alpha_zero;

      // check diagonal entries
      if (i == j)
      {
        REQUIRE(Approx(result).margin(toterr) == asymp);
        std::cout << std::scientific;
        std::cout << result << std::endl;
        std::cout << asymp << std::endl;
      }
      else
      {
        REQUIRE(result == 0); // off-diagonals are zero
      };

    };
  };
};
