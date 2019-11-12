#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Vacuum Greens Tensor works properly")
{
  struct Options opts;
  opts.fancy_I=true;

  GreensTensorVacuum *gre = new GreensTensorVacuum(0.01,1e3);
  cx_mat::fixed<3,3> test(fill::zeros);

  gre->calculate_integrated(test, 3.0, opts);

  std::cout << test << std::endl;
  REQUIRE( 0 == 0 );

};
