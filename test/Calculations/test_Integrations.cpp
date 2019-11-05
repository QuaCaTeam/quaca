#include "catch.hpp"
#include "Quaca.h"
#include <iostream>
#include <math.h>

double f( double x, void *p)
{
  return 1E0/(x*x+1E0);
};

double f2( double x, void *p)
{
  return 1E0/sqrt(1E0-x);
};

TEST_CASE("Integration test")
{
  SECTION("CQUAD yields the demanded accuracy")
  {
    double testcquad=cquad(f,-1E10,1E10,1E-10,0);
    REQUIRE( testcquad == Approx(M_PI).epsilon(1E-10));
  };
  SECTION("QAGS yields the demanded accuracy")
  {
    double testqags=qags(f2,0E0,1E0,1E-10,0);
    REQUIRE( testqags == Approx(2E0).epsilon(1E-10));
  };
  SECTION("QAGIU yields the demanded accuracy")
  {
    double testqagiu=qagiu(f,0,1E-10,0);
    REQUIRE( testqagiu == Approx(M_PI/2E0).epsilon(1E-10));
  };

};
