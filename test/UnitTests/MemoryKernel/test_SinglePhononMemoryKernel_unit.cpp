#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("SinglePhonon memory kernel constructors work as expected",
          "[SinglePhononMemoryKernel]") {

  SECTION("Direct constructor") {
    SinglePhononMemoryKernel memorykernel(1e-1,1.e-5,4.34,1e-5);
    REQUIRE(memorykernel.get_gamma() == 1e-1);
    REQUIRE(memorykernel.get_gamma_phon() == 1e-5);
    REQUIRE(memorykernel.get_omega_phon() == 4.34);
    REQUIRE(memorykernel.get_coupling() == 1e-5);
  }

  SECTION("json file constructor") {
    SinglePhononMemoryKernel memorykernel("../data/test_files/SinglePhononMemoryKernel.json");
    REQUIRE(memorykernel.get_gamma() == 1e-1);
    REQUIRE(memorykernel.get_gamma_phon() == 1e-5);
    REQUIRE(memorykernel.get_omega_phon() == 4.34);
    REQUIRE(memorykernel.get_coupling() == 1e-5);
  }
}

TEST_CASE("SinglePhonon memory kernel obeys crossing relation",
          "[SinglePhononMemoryKernel]") {
  SinglePhononMemoryKernel mk("../data/test_files/SinglePhononMemoryKernel.json");
  REQUIRE(mk.calculate(1.0) == std::conj(mk.calculate(-1.0)));
}
