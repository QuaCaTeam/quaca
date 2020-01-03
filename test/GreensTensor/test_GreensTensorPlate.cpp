#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Vacuum Greens Tensor works properly")
{
  struct Options_GreensTensor opts;
  opts.fancy_I=true;
  opts.omega = 3.0;

  GreensTensorVacuum Greens(0.01,1e3);
  opts.class_pt = &Greens;

  cx_mat::fixed<3,3> test(fill::zeros);

  Greens.integrate_k_1d(test, opts);

  double analytic_result_prefactor = 2./3.*pow(opts.omega,3)/pow(1-pow(Greens.get_v(),2),2);

  std::cout << test << std::endl;
  REQUIRE( Approx(test(0,0).real()).epsilon(1E-4) == analytic_result_prefactor );
  REQUIRE( Approx(test(1,1).real()).epsilon(1E-4) == analytic_result_prefactor*(1+pow(Greens.get_v(),2))/(1-pow(Greens.get_v(),2)));
  REQUIRE( Approx(test(2,2).real()).epsilon(1E-4) == analytic_result_prefactor*(1+pow(Greens.get_v(),2))/(1-pow(Greens.get_v(),2)));

};
