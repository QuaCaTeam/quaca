#include "Quaca.h"
#include "catch.hpp"
#include <iostream>

TEST_CASE("Options constructors work as expected", "[Options]") {
  SECTION("Direct constructor") {
    std::string file = "../README.json";
    unsigned int threads = 3;
    Options opts(file, threads);
    REQUIRE(opts.get_parameter_file() == file);
    REQUIRE(opts.get_output_file() == "../README.csv");
    REQUIRE(opts.get_num_threads() == 3);
  }

}
