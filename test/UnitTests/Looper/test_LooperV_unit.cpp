#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LooperV constructors work as expected", "[LooperV]") {
  SECTION("Direct constructor") {
    double start = 1e-3;
    double end = 1e-1;
    int number_of_steps = 10;
    std::string scale = "log";

    LooperV looper(start, end, number_of_steps, scale);

    REQUIRE(looper.get_steps_total() == number_of_steps);
    REQUIRE(looper.get_step(0) == start);
    REQUIRE(looper.get_step(number_of_steps - 1) == Approx(end));
  }

  SECTION("json file constructor") {
    LooperV looper("../data/test_files/LooperV.json");

    REQUIRE(looper.get_steps_total() == 20);
    REQUIRE(looper.get_step(0) == 10.2);
    REQUIRE(looper.get_step(19) == 35.5);
  }
}

TEST_CASE("LooperV Steps are calculated correctly", "[LooperV]") {
  SECTION("Steps for linear scale") {
    double start = 0;
    double end = 3;
    int number_of_steps = 4;
    std::string scale = "linear";

    LooperV looper(start, end, number_of_steps, scale);

    REQUIRE(looper.get_step(0) == Approx(start));
    REQUIRE(looper.get_step(1) == Approx(1));
    REQUIRE(looper.get_step(2) == Approx(2));
    REQUIRE(looper.get_step(3) == Approx(end));
  }

  SECTION("Steps for logarithmic scale") {
    double start = 1e-4;
    double end = 1e-1;
    int number_of_steps = 4;
    std::string scale = "log";

    LooperV looper(start, end, number_of_steps, scale);

    REQUIRE(looper.get_step(0) == Approx(start));
    REQUIRE(looper.get_step(1) == Approx(1e-3));
    REQUIRE(looper.get_step(2) == Approx(1e-2));
    REQUIRE(looper.get_step(3) == Approx(end));
  }
}
