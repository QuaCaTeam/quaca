#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LooperZa constructors work as expected", "[LooperZa]") {
  SECTION("Direct constructor") {
    double start = 0;
    double end = 1;
    int number_of_steps = 10;
    std::string scale = "log";

    LooperZa looper(start, end, number_of_steps, scale);

    REQUIRE(looper.get_steps_total() == number_of_steps);
    REQUIRE(looper.get_step(0) == start);
    REQUIRE(looper.get_step(number_of_steps - 1) == end);
  };

  SECTION("ini file constructor") {
    LooperZa looper("../data/test_files/LooperZa.ini");

    REQUIRE(looper.get_steps_total() == 20);
    REQUIRE(looper.get_step(0) == 10.2);
    REQUIRE(looper.get_step(19) == 35.5);
  };
};

TEST_CASE("LooperZa Steps are calculated correctly", "[LooperZa]") {
  SECTION("Steps for linear scale") {
    double start = 0;
    double end = 3;
    int number_of_steps = 4;
    std::string scale = "linear";

    LooperZa looper(start, end, number_of_steps, scale);

    REQUIRE(looper.get_step(0) == Approx(start));
    REQUIRE(looper.get_step(1) == Approx(1));
    REQUIRE(looper.get_step(2) == Approx(2));
    REQUIRE(looper.get_step(3) == Approx(end));
  };

  SECTION("Steps for logarithmic scale") {
    double start = 1e-4;
    double end = 1e-1;
    int number_of_steps = 4;
    std::string scale = "log";

    LooperZa looper(start, end, number_of_steps, scale);

    REQUIRE(looper.get_step(0) == Approx(start));
    REQUIRE(looper.get_step(1) == Approx(1e-3));
    REQUIRE(looper.get_step(2) == Approx(1e-2));
    REQUIRE(looper.get_step(3) == Approx(end));
  };
};
