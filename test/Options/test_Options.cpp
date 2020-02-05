#include "Quaca.h"
#include "catch.hpp"
#include <iostream>

TEST_CASE("Options constructors work", "[Options]") {
  SECTION("direct construction") {
    std::string file = "../README.ini";
    Options opts(file);
    REQUIRE(opts.get_parameter_file() == file);
    REQUIRE(opts.get_output_file() == "../README.csv");
  };

  SECTION("construction with command line arguments") {
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
