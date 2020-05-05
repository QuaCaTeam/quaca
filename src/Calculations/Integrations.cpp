#include "Integrations.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

// wrapper to cquad routine
double cquad(const std::function<double(double)> &f, double a, double b,
             double relerr, double epsabs) {

  // cast given lambda to static (!) gsl_function
  gsl_function_pp<decltype(f)> Fp(f);
  auto *F = static_cast<gsl_function *>(&Fp);

  double res;
  /* Initialize the workspace. */
  gsl_integration_cquad_workspace *ws =
      gsl_integration_cquad_workspace_alloc(100);
  if (ws == nullptr) {
    printf("call to gsl_integration_cquad_workspace_alloc failed.\n");
    abort();
  }

  /* Call the integrator. */
  /* set nevals and abserr pointer to nullptr, we are only interested in result
   */
  int success = gsl_integration_cquad(F, a, b, epsabs, relerr, ws, &res,
                                      nullptr, nullptr);
  if (success != 0) {
    printf("cquad error: %s\n", gsl_strerror(success));
    abort();
  }

  /* Free the workspace. */
  gsl_integration_cquad_workspace_free(ws);

  return res;
}

double qags(const std::function<double(double)> &f, double a, double b,
            double relerr, double epsabs) {

  // cast given lambda to static (!) gsl_function
  gsl_function_pp<decltype(f)> Fp(f);
  auto *F = static_cast<gsl_function *>(&Fp);

  double res;
  double abserr;

  /* Initialize the workspace. */
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);
  if (ws == nullptr) {
    printf("call to gsl_integration_workspace_alloc failed.\n");
    abort();
  }

  /* Call the integrator. */
  /* set nevals and abserr pointer to nullptr, we are only interested in result
   */
  int success =
      gsl_integration_qags(F, a, b, epsabs, relerr, 1000, ws, &res, &abserr);
  if (success != 0) {
    printf("qags error: %s\n", gsl_strerror(success));
    abort();
  }

  /* Free the workspace. */
  gsl_integration_workspace_free(ws);

  return res;
}

// wrapper to qagiu routine
double qagiu(const std::function<double(double)> &f, double a, double relerr,
             double epsabs) {

  // cast given lambda to static (!) gsl_function
  gsl_function_pp<decltype(f)> Fp(f);
  auto *F = static_cast<gsl_function *>(&Fp);

  double res;
  double abserr;

  /* Initialize the workspace. */
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);
  if (ws == nullptr) {
    printf("call to gsl_integration_workspace_alloc failed.\n");
    abort();
  }

  /* Call the integrator. */
  /* set nevals and abserr pointer to nullptr, we are only interested in result
   */
  int success =
      gsl_integration_qagiu(F, a, epsabs, relerr, 1000, ws, &res, &abserr);
  if (success != 0) {
    printf("qagiu error: %s\n", gsl_strerror(success));
    abort();
  }

  return res;
}
