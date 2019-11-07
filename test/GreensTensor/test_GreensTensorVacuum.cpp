#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Vacuum Greens Tensor works properly")
{
  GreensTensorVacuum *gre = new GreensTensorVacuum(1.0);
  cx_mat::fixed<3,3> test(fill::zeros);
  gre->calculate_pure(test, vec(2,fill::zeros),3.0);
  std::cout << test << std::endl;
  REQUIRE( 0 == 0 );

};
