#include "Quaca.h"
#include "catch.hpp"
#include <iostream>

TEST_CASE("Options constructors work as expected", "[Options]") {
  SECTION("Direct constructor") {
    std::string file = "../README.ini";
    Options opts(file);
    REQUIRE(opts.get_parameter_file() == file);
    REQUIRE(opts.get_output_file() == "../README.csv");
  };

  SECTION("Command line argument constructor") {
    // simulate command line arguments
    int argc = 3;
    char *argv[argc];
    argv[0] = (char *)"./../bin/Quaca";
    argv[1] = (char *)"--file";
    argv[2] = (char *)"../data/NonExistantFile.ini";

    Options opts(argc, argv);

    REQUIRE(opts.get_parameter_file() == "../data/NonExistantFile.ini");
    REQUIRE(opts.get_output_file() == "../data/NonExistantFile.csv");
  };
};
