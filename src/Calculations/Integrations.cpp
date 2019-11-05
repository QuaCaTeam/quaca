#include "Integrations.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>


double cquad(double my_f(double, void *), double a , double b, double relerr, double epsabs) {
    double res;

    /* Prepare the function. */
    gsl_function f;
    f.function = my_f;

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
