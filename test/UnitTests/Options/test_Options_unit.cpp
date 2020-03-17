#include "Quaca.h"
#include "catch.hpp"
#include <iostream>

TEST_CASE("Options constructors work as expected", "[Options]") {
  SECTION("Direct constructor") {
    std::string file = "../README.json";
    unsigned int threads = 3;
    Options opts(file, 3);
    REQUIRE(opts.get_parameter_file() == file);
    REQUIRE(opts.get_output_file() == "../README.csv");
    REQUIRE(opts.get_num_threads() == 3);
  };

  SECTION("Command line argument constructor") {
    // simulate command line arguments
    int argc = 5;
    char *argv[argc];
    argv[0] = (char *)"./../bin/Quaca";
    argv[1] = (char *)"--file";
    argv[2] = (char *)"../data/NonExistantFile.json";
    argv[3] = (char *)"--threads";
    argv[4] = (char *)"3";

    Options opts(argc, argv);

    REQUIRE(opts.get_parameter_file() == "../data/NonExistantFile.json");
    REQUIRE(opts.get_output_file() == "../data/NonExistantFile.csv");
    REQUIRE(opts.get_num_threads() == 3);
  };
};
