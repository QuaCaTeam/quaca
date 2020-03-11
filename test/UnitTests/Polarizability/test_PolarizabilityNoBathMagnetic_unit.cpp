//
// Created by hermasim on 11/03/2020.
//

#include <armadillo>
#include <complex>

#include "Quaca.h"
#include "catch.hpp"

/*
TEST_CASE("Ensure that PolarizabilityNoBathMagnetic fulfills the crossing relation"
          "[PolarizabilityNoBath]") {
  //define Green's tesnsor
  double omega_p = 9;
  double gamma = 0.1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic greens(v, za, 0.1, &refl, delta_cut, rel_err);


  // define polarizability
  auto omega_a = GENERATE(take(1, random(1e-1, 1e1)));
  auto alpha_zero = GENERATE(take(1, random(0.0, 0.1)));
  PolarizabilityNoBathMagnetic pol(omega_a, alpha_zero, &greens);

  auto omega = GENERATE(take(3, random(0.,1e2)));

  Options_Polarizability opts;
  opts.class_pt = &pol;

  //Matrices to store the results
  cx_mat::fixed<3,3> LHS(fill::zeros);
  cx_mat::fixed<3,3> RHS(fill::zeros);

  opts.omega = omega;
  pol.calculate_tensor(LHS, opts);
  opts.omega = - omega;
  pol.calculate_tensor(RHS, opts);
  REQUIRE(approx_equal(LHS, conj(RHS), "abs", 1e-12));
};
*/
/*
TEST_CASE("Ensure that there are non negligible contributions in the magnetic "
          "part of the polarizability"
          "[PolarizabilityNoBath]") {
  //define Green's tesnsor
  double omega_p = 9;
  double gamma = 0.1;
  double v = 1e-2;
  double za = 0.1;
  double delta_cut = 30;
  vec::fixed<2> rel_err = {1E-8, 1E-6};
  PermittivityDrude perm(omega_p, gamma);
  ReflectionCoefficientsLocBulk refl(&perm);
  GreensTensorPlateMagnetic greens_magnetic(v, za, 0.1, &refl, delta_cut, rel_err);
  GreensTensorPlate greens(v, za, 0.1, &refl, delta_cut, rel_err);


  // define polarizability
  auto omega_a = GENERATE(take(1, random(1e-1, 1e1)));
  auto alpha_zero = GENERATE(take(1, random(0.0, 0.1)));
  PolarizabilityNoBathMagnetic pol_magnetic(omega_a, alpha_zero, &greens_magnetic);
  PolarizabilityNoBath pol(omega_a, alpha_zero, &greens);

  auto omega = GENERATE(take(3, random(0.,1e2)));

  Options_Polarizability opts;

  //Matrices to store the results
  cx_mat::fixed<3,3> LHS(fill::zeros);
  cx_mat::fixed<3,3> RHS(fill::zeros);

  opts.omega = omega;
  //compute Polarizability with electric Green's tensor only
  opts.class_pt = &pol;
  pol.calculate_tensor(LHS, opts);

  //compute Polarizability with magnetic Green's tensor as well
  opts.class_pt = &pol_magnetic;
  pol_magnetic.calculate_tensor(RHS, opts);

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      REQUIRE(LHS(i,j) != RHS(i,j));
    }
  }
};
 */
