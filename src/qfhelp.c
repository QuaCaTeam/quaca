/* --------- */
/* LIBRARIES */
/* --------- */

#include "h/qfhelp.h"


/* --------- */
/* FUNCTIONS */
/* --------- */

// matrix multiplication
void multiply(double complex mat1[Ndim][Ndim], double complex mat2[Ndim][Ndim], double complex res[Ndim][Ndim]) {
    int i, j, k;
    for (i = 0; i < Ndim; i++) {
        for (j = 0; j < Ndim; j++) {
            res[i][j] = 0.;
            for (k = 0; k < Ndim; k++) {
                res[i][j] += mat1[i][k]*mat2[k][j]; // add entries
            }
        }
    }
};

// complex transpose
void transpose(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim]) {
    int i, j;
    for (i = 0; i < Ndim; i++) {
        for (j = 0; j < Ndim; j++) {
            res[j][i] = mat[i][j]; // transpose
        }
    }
};

//  transpose and complex conjugate
void dagger(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim]) {
    int i, j;
    for (i = 0; i < Ndim; i++) {
        for (j = 0; j < Ndim; j++) {
            res[j][i] = conj(mat[i][j]); // transpose and complex conjugate
        }
    }
};

// fancy I or R
void fancy(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim], int mode) {
    int i, j;
    for (i = 0; i < Ndim; i++) {
        for (j = 0; j < Ndim; j++) {
            res[i][j] = -0.5E0*I*( mat[i][j] + mode*conj(mat[j][i]) );
        }
    }
};

// trace of matrix
const double complex tr(double complex mat[Ndim][Ndim]) {
    int i;
    double complex res;
    res = 0.;
    for (i = 0; i < Ndim; i++) {
        res += mat[i][i]; // add diagonal entries
    }
    return res;
};

double cquad(double my_f(), void * p, double a , double b, double relerr, double epsabs) {
    gsl_function f;
    gsl_integration_cquad_workspace *ws = NULL;
    struct parameters * params = (struct parameters *)p;
    double res, abserr;
    size_t neval;

    /* Prepare the function. */
    f.function = my_f;
    f.params = params;

    /* Initialize the workspace. */
    if ( ( ws = gsl_integration_cquad_workspace_alloc( 200 ) ) == NULL ) {
        printf( "call to gsl_integration_cquad_workspace_alloc failed.\n" );
        abort();
    }

    /* Call the integrator. */
    if ( gsl_integration_cquad( &f, a , b , epsabs , relerr , ws , &res , &abserr , &neval ) != 0 ) {
        printf( "call to gsl_integration_cquad failed.\n" );
        abort();
    }

    /* Free the workspace. */
    gsl_integration_cquad_workspace_free( ws );

    return res;
};

double qags(double my_f(), void * p, double a , double b, double relerr, double epsabs) {
    gsl_function f;
    gsl_integration_workspace *ws = NULL;
    struct parameters * params = (struct parameters *)p;
    double res, abserr;
    size_t limit = 1000;

    /* Prepare the function. */
    f.function = my_f;
    f.params = params;

    /* Initialize the workspace. */
    if ( ( ws = gsl_integration_workspace_alloc( 1000 ) ) == NULL ) {
        printf( "call to gsl_integration_cquad_workspace_alloc failed.\n" );
        abort();
    }

    /* Call the integrator. */
    if ( gsl_integration_qags( &f, a , b , epsabs , relerr , limit, ws , &res , &abserr ) != 0 ) {
        printf( "call to gsl_integration_cquad failed.\n" );
        abort();
    }

    /* Free the workspace. */
    gsl_integration_workspace_free( ws );

    return res;
};

double integinf(double my_f(), void * p, double a, double relerr, double epsabs) {
    gsl_function f;
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);
    struct parameters * params = (struct parameters *)p;
    double res, abserr;

    /* Prepare the function. */
    f.function = my_f;
    f.params = params;


    /* Call the integrator. */
    if ( gsl_integration_qagiu( &f, a, epsabs, relerr, 1000, work_ptr, &res, &abserr) != 0 ) {
        printf( "call to gsl_integration_qagiu failed.\n" );
        abort();
     }

    /* Free the workspace. */
    gsl_integration_workspace_free( work_ptr );

    return res;
};
