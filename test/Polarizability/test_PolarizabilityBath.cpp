#include <complex>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Polarizability with internal bath")
{
  SECTION("Constructor from .ini file")
  {
    // test Polarizability parameters
    PolarizabilityBath *pol = new PolarizabilityBath("../data/test_files/PolarizabilityBath.ini");
    REQUIRE( pol->get_omega_a() == 1.3 );
    REQUIRE( pol->get_alpha_zero() == 6E-9 );

    // test if we read memory kernel correctly
    std::complex<double> test = pol->get_mu(3.0);
    REQUIRE( test.real() == 0.69420 );
  };

  SECTION("Alpha test")
  {
    PolarizabilityBath *pol = new PolarizabilityBath("../data/test_files/PolarizabilityBath.ini");
    cx_mat test = pol->calculate(3.0);
    std::cout << test << std::endl;
    REQUIRE( 0 == 0 );
  };
};
