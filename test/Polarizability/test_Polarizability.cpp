#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("MuTest")
{
  Polarizability *pol = new Polarizability();
SECTION( "mu(w) has an even symmetry due to causality" ){
    REQUIRE( pol->muR(-1E0) == pol->muR(1E0));
};
SECTION( "mu(w) is 2" ){
    REQUIRE( pol->muR(-1E0) == 2);
};
};
