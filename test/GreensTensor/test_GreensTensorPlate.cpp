#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Plate Greens Tensor works properly")
{
  struct Options_GreensTensor opts;
  opts.fancy_I=true;
  opts.omega = 3.0;
  opts.kvec = { 10.0 , NAN };
  opts.indices = { 0 , 0};
  GreensTensorPlate Greens(0.01 , 1e-1, 1e3, "../data/test_files/PermittivityDrude.ini");
  opts.class_pt = &Greens;

//  std::cout << Greens.integrand_k_2d( 1., &opts) << std::endl;
  std::cout << Greens.integrand_k_1d( 1., &opts) << std::endl;
};
