#include "Quaca.h"
#include "catch.hpp"
#include <cmath>

TEST_CASE("Integration routines return right results", "[Integrations]") {
  SECTION("CQUAD yields the demanded accuracy") {
    auto f = [=](double x) -> double {return 1./(x * x + 1.0);};
    double testcquad = cquad(f, -1E10, 1E10, 1E-10, 0);
    REQUIRE(testcquad == Approx(M_PI).epsilon(1E-10));
  }

  SECTION("QAGS yields the demanded accuracy") {
    auto g = [=](double x) -> double {return 1./sqrt(1.0 - x);};
    double testqags = qags(g, 0.0, 1.0, 1E-10, 0);
    REQUIRE(testqags == Approx(2.0).epsilon(1E-10));
  }

  SECTION("QAGIU yields the demanded accuracy") {
    auto f = [=](double x) -> double {return 1./(x * x + 1.0);};
    double testqagiu = qagiu(f, 0, 1E-10, 0);
    REQUIRE(testqagiu == Approx(M_PI / 2.0).epsilon(1E-10));
  }
}
