#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("Lorentz memory kernel constructors work as expected",
          "[LorentzMemoryKernel]") {

  SECTION("Direct constructor") {
    LorentzMemoryKernel memorykernel(1.e-5,4.34,14.28,1.035);
    REQUIRE(memorykernel.get_gamma() == 1e-5);
    REQUIRE(memorykernel.get_omega_0() == 4.34);
    REQUIRE(memorykernel.get_omega_p() == 14.28);
    REQUIRE(memorykernel.get_eps_inf() == 1.035);
  }

  SECTION("json file constructor") {
    LorentzMemoryKernel memorykernel("../data/test_files/LorentzMemoryKernel.json");
    REQUIRE(memorykernel.get_gamma() == 1e-5);
    REQUIRE(memorykernel.get_omega_0() == 4.34);
    REQUIRE(memorykernel.get_omega_p() == 14.28);
    REQUIRE(memorykernel.get_eps_inf() == 1.035);
  }
}

TEST_CASE("Lorentz memory kernel obeys crossing relation",
          "[LorentzMemoryKernel]") {
  LorentzMemoryKernel mk("../data/test_files/LorentzMemoryKernel.json");
  REQUIRE(mk.calculate(1.0) == std::conj(mk.calculate(-1.0)));
}
