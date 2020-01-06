#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Plate Greens Tensor works properly")
{
  struct Options_GreensTensor opts;
  opts.fancy_I=true;
  opts.omega = 3.0;
  opts.kvev = { 10.0 , NAN };

  GreensTensorPlate Greens(0.01 , 1e-1, 1e3, "../data/test_files/PermittivityDrude.ini");
  opts.class_pt = &Greens;

  std::cout << Greens.integrate_k_2d( 1, opts) << std::endl;
};
