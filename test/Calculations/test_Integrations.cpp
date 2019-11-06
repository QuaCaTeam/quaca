#include "catch.hpp"
#include "Quaca.h"
#include <iostream>
#include <math.h>

// f(x) = 1/(x^2 +1)
double f( double x, void *p)
{
  return 1.0/(x*x+1.0);
};

// g(x) = 1/(1-x)^(1/2);
double g( double x, void *p)
{
  return 1.0/sqrt(1.0-x);
};

TEST_CASE("Integration routines return right results", "[Integrations]")
{
  SECTION("CQUAD yields the demanded accuracy")
  {
    double testcquad = cquad(f, -1E10, 1E10, 1E-10, 0);
    REQUIRE( testcquad == Approx(M_PI).epsilon(1E-10) );
  };

  SECTION("QAGS yields the demanded accuracy")
  {
    double testqags = qags(g, 0.0, 1.0, 1E-10, 0);
    REQUIRE( testqags == Approx(2.0).epsilon(1E-10) );
  };

  SECTION("QAGIU yields the demanded accuracy")
  {
    double testqagiu =qagiu(f, 0, 1E-10, 0);
    REQUIRE( testqagiu == Approx(M_PI/2.0).epsilon(1E-10) );
  };

};
