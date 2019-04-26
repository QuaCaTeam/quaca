#include <check.h>
#include "../src/plate.h"
#include "../src/qfhelp.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <string.h>

/* General Tests */
START_TEST(test_trivial) {
  struct parameters inputparams;
  inputparams.gamMu = 3.0;
  double complex test = mu(3.0, &inputparams);
  ck_assert_double_eq(3.0, creal(test));
} END_TEST


START_TEST(test_trivial_two) {
  struct parameters inputparams;
  inputparams.gamMu = 3.0;
  double complex test = mu(3.0, &inputparams);
  ck_assert_double_eq(3.001, creal(test));
} END_TEST

Suite *trivial_suite(void) {
  Suite *s;
  TCase *tc_core;
  TCase *tc_second;

  s = suite_create("Trivial Suite");

  tc_core = tcase_create("First");
  tcase_add_test(tc_core, test_trivial_two);
  suite_add_tcase(s, tc_core);

  tc_second = tcase_create("Second");
  tcase_add_test(tc_second, test_trivial);
  suite_add_tcase(s, tc_second);

  return s;
}

/* reflection coefficients */
START_TEST(ref_coeff_test) {
  struct parameters inputparams;
  inputparams.einf = 1.0;
  inputparams.wp1 = 9;
  inputparams.g1 = 30E-3;

  double complex ref[2], ref2[2];
  refl(ref, 3, 4, &inputparams);
  refl(ref2, -3, 4, &inputparams);
  ck_assert_double_eq(cimag(ref[0]), -cimag(ref2[0]));
  ck_assert_double_eq(cimag(ref[1]), -cimag(ref2[1]));
} END_TEST

Suite *ref_suite(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Reflection Coefficients");
  tc_core = tcase_create("Number check");

  tcase_add_test(tc_core, ref_coeff_test);
  suite_add_tcase(s, tc_core);

  return s;
}

/* polarizability*/
START_TEST(polar_test) {
  struct parameters inputparams;
  inputparams.v = 5E-4;
  inputparams.za = 5E-9/(6.58E-16*2.99E8);
  inputparams.muquest = 1;
  inputparams.gamMu = 1E-1;
  inputparams.wa = 1.3;
  inputparams.a0 = 6E-9;

  inputparams.eps0 = 1/(4*PI);

  inputparams.einf = 1;
  inputparams.wp1 = 9;
  inputparams.g1 = 35E-3;
  inputparams.beta = 1/(1E-6*8.6173303E-5);

  inputparams.kcut = 30;
  inputparams.relerr = 1E-2;
  inputparams.recerr = 1E-2;
  inputparams.abserr = 1E-200;

  double complex alp[3][3], alp2[3][3];
  alpha(alp, 3, &inputparams);
  alpha(alp2, -3, &inputparams);

  ck_assert_double_eq_tol(cimag(alp[0][0]), -cimag(alp2[0][0]), 1E-13);
} END_TEST

Suite *polar_suite(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Polarizability");
  tc_core = tcase_create("Number check");

  tcase_add_test(tc_core, polar_test);
  suite_add_tcase(s, tc_core);

  return s;
}


int main(void) {

    int no_failed = 0;
   
    // create Suites 
    Suite *trivial_suite (void);
    Suite *other_suite (void);
    SRunner *runner;

    // add the suites
    runner = srunner_create(trivial_suite ());
    srunner_add_suite(runner, ref_suite ());
    srunner_add_suite(runner, polar_suite ());

    // run all tests
    srunner_set_log(runner, "test/test.log");
    srunner_run_all(runner, CK_NORMAL);
    no_failed = srunner_ntests_failed(runner);
    srunner_free(runner);
    return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
