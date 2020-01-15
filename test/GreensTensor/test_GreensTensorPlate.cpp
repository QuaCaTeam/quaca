#include <complex>
#include <armadillo>

#include "catch.hpp"
#include "Quaca.h"

TEST_CASE("Plate Greens Tensor works properly")
{
  struct Options_GreensTensor opts;
  opts.fancy_R=true;
  opts.omega = 3.0;
  opts.kvec = { 10.0 , NAN };
  opts.indices = { 0 , 0};
  GreensTensorPlate Greens(0.01 , 1e-1, 1e3, "../data/test_files/GreensTensorPlate.ini");
  opts.class_pt = &Greens;

  cx_mat::fixed<3,3> test(fill::zeros);

  Greens.integrate_k_1d(test, opts);
  std::cout << test << std::endl;
//  std::cout << Greens.integrand_k_2d( 1., &opts) << std::endl;
  std::cout << Greens.integrand_k_1d( 1., &opts) << std::endl;
};

TEST_CASE("Vacuum Greens Tensor works properly")
{

  GreensTensorVacuum Greens(0.01,1e3);
  /*!
   * Some basic relations any Green's tensor should fulfill which can
   * be found in docs under: Relations_and_tests.pdf
   */
  vec::fixed<2> kvec = ones<vec>(2);
  double omega = 5;
  cx_mat::fixed<3,3> Greens_lhs(fill::zeros);
  cx_mat::fixed<3,3> Greens_rhs(fill::zeros);
  /*!
   * Crossing relation in frequency domain see eq. [1]
   */
  Greens.calculate_tensor(Greens_lhs, -kvec, -omega);
  Greens.calculate_tensor(Greens_rhs, kvec, omega);
  REQUIRE(approx_equal(Greens_lhs, trans(conj(Greens_rhs)),"absdiff", 10E-5));

  /*!
   * Reciprocity, see eq. [6]
   */
  Greens.calculate_tensor(Greens_lhs, -kvec, omega);
  Greens.calculate_tensor(Greens_rhs, kvec, omega);
  REQUIRE(approx_equal(Greens_lhs, trans(Greens_rhs),"absdiff", 10E-5));

  /*!
   * Reality, see eq. [7], due to the divergence of the real part only the imaginary
   * part is calculated in any case
   */
  Greens.calculate_tensor(Greens_lhs, kvec, omega);
  Greens.calculate_tensor(Greens_rhs, kvec, -omega);
  REQUIRE(approx_equal(Greens_lhs, Greens_rhs,"absdiff", 10E-5));



  /*!
   * Test the integration routine of the vacuums Green's tensor, the analytical results
   * can be found in the docs under: VacuumGreen.pdf, see eq. [21]
   */
  struct Options_GreensTensor opts;
  opts.fancy_I=true;
  opts.omega = 3.0;
  opts.class_pt = &Greens;

  cx_mat::fixed<3,3> test(fill::zeros);

  Greens.integrate_k_1d(test,opts);

  double analytic_result_prefactor = 2./3.*pow(opts.omega,3)/pow(1-pow(Greens.get_v(),2),2);

  std::cout << test << std::endl;
  REQUIRE( Approx(test(0,0).real()).epsilon(1E-4) == analytic_result_prefactor );
  REQUIRE( Approx(test(1,1).real()).epsilon(1E-4) == analytic_result_prefactor*(1+pow(Greens.get_v(),2))/(1-pow(Greens.get_v(),2)));
  REQUIRE( Approx(test(2,2).real()).epsilon(1E-4) == analytic_result_prefactor*(1+pow(Greens.get_v(),2))/(1-pow(Greens.get_v(),2)));

};
