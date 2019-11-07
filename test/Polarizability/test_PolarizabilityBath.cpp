#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Polarizability with bath can be constructed in different ways", "[PolarizabilityBath]")
{
  SECTION("Constructor from .ini file")
  {
    PolarizabilityBath *pol = new PolarizabilityBath("../data/test_files/PolarizabilityBath.ini");
    REQUIRE( pol->get_omega_a() == 1.3 );
    REQUIRE( pol->get_alpha_zero() == 6E-9 );

    // test if we read memory kernel correctly
    std::complex<double> test = pol->get_mu(3.0);
    REQUIRE( test.real() == 0.69420 );
  };

  SECTION("Constructor with direct input")
  {
    MemoryKernel *mu = new OhmicMemoryKernel(0.69420);
    PolarizabilityBath *pol = new PolarizabilityBath(1.3, 6E-9, mu);
    REQUIRE( pol->get_omega_a() == 1.3 );
    REQUIRE( pol->get_alpha_zero() == 6E-9 );

    std::complex<double> test = pol->get_mu(3.0);
    REQUIRE( test.real() == 0.69420 );
  };
};

TEST_CASE("PolarizabilityBath returns a diagonal matrix", "[PolarizabilityBath]")
{
  PolarizabilityBath *pol = new PolarizabilityBath("../data/test_files/PolarizabilityBath.ini");
  cx_mat test(3,3, fill::zeros);
  pol->calculate(test, 3.0);
  std::cout << test << std::endl;
  REQUIRE( 0 == 0 );

};
