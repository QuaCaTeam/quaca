#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Construction of Green's tensor works properly", "[GreensTensorPlate]")
{
  SECTION("Construction with direct input")
  {
  double omega   = 1;
  double omega_p = 9;
  double gamma   = 0.1;
  double v = 1E-5;
  double za = 0.1;
  double beta = 1E4;
  double delta_cut = 20;

  PermittivityDrude perm(gamma, omega_p);
  GreensTensorPlate Greens(v , za, beta, &perm, delta_cut);
  REQUIRE(Greens.get_za() == za);
  REQUIRE(Greens.get_delta_cut() == delta_cut);
  };

  SECTION("Construction from .ini file")
  {
  GreensTensorPlate Greens("../data/test_files/GreensTensorPlate.ini");
  REQUIRE(Greens.get_za() == 0.1);
  REQUIRE(Greens.get_delta_cut() == 10);
  };
};
TEST_CASE("The operations calculate_tensor and the integrand_k_2d coincide", "[GreensTensorPlate]")
{

};

TEST_CASE("Scattered Green's tensor works properly", "[GreensTensorPlate]")
{

};

//  struct Options_GreensTensor opts;
//  opts.fancy_R=true;
//  opts.omega = 3.0;
//  opts.kvec = { 10.0 , NAN };
//  opts.indices = { 0 , 0};
//  opts.class_pt = &Greens;
//  cx_mat::fixed<3,3> test(fill::zeros);
//  Greens.integrate_k_1d(test, opts);
//  std::cout << test << std::endl;
////  std::cout << Greens.integrand_k_2d( 1., &opts) << std::endl;
//  std::cout << Greens.integrand_k_1d( 1., &opts) << std::endl;
