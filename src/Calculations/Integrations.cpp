#include "Integrations.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

// wrapper to cquad routine
double cquad(double my_f(double, void *), void* params, double a , double b, double relerr, double epsabs)
{
  double res;

  /* Prepare the function. */
  gsl_function f;
  f.function = my_f;
  f.params = params;

  /* Initialize the workspace. */
  gsl_integration_cquad_workspace *ws = gsl_integration_cquad_workspace_alloc(100);
  if ( ws == NULL ) {
    printf( "call to gsl_integration_cquad_workspace_alloc failed.\n" );
    abort();
  }

  /* Call the integrator. */
  /* set nevals and abserr pointer to NULL, we are only interested in result */
  int success = gsl_integration_cquad( &f, a , b , epsabs , relerr , ws , &res , NULL , NULL );
  if ( success != 0 ) {
    printf( "call to gsl_integration_cquad failed.\n" );
    abort();
  }

  /* Free the workspace. */
  gsl_integration_cquad_workspace_free( ws );

  return res;
};

// wrapper to qags routine
double qags(double my_f(double, void *), void* params, double a , double b, double relerr, double epsabs)
{
  double res;
  double abserr;

  /* Prepare the function. */
  gsl_function f;
  f.function = my_f;
  f.params = params;

  /* Initialize the workspace. */
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);
  if ( ws == NULL ) {
    printf( "call to gsl_integration_workspace_alloc failed.\n" );
    abort();
  }

  /* Call the integrator. */
  /* set nevals and abserr pointer to NULL, we are only interested in result */
  int success = gsl_integration_qags( &f, a , b , epsabs , relerr , 1000 ,ws , &res , &abserr );
  if ( success != 0 ) {
    printf( "call to gsl_integration_qags failed.\n" );
    abort();
  }

  /* Free the workspace. */
  gsl_integration_workspace_free( ws );

  return res;

};

// wrapper to qagiu routine
double qagiu(double my_f(double, void *), void* params, double a , double relerr, double epsabs)
{
  double res;
  double abserr;

  /* Prepare the function. */
  gsl_function f;
  f.function = my_f;
  f.params = params;
  
  /* Initialize the workspace. */
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);
  if ( ws == NULL ) {
    printf( "call to gsl_integration_workspace_alloc failed.\n" );
    abort();
  }

  /* Call the integrator. */
  /* set nevals and abserr pointer to NULL, we are only interested in result */
  int success = gsl_integration_qagiu( &f, a , epsabs , relerr , 1000 ,ws , &res , &abserr );
  if ( success != 0 ) {
    printf( "call to gsl_integration_qagiu failed.\n" );
    abort();
  }

  return res;

  /* Free the workspace. */
  gsl_integration_workspace_free( ws );

};
