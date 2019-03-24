/* --------- */
/* LIBRARIES */
/* --------- */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <string.h>


/* ---------- */
/* PARAMETERS */
/* ---------- */

/* mathematical constants */
#define PI 3.14159265358979323846
#define Ndim 3

/* physical constants */
#define kB 8.6173303E-5
#define hbar 6.582119514E-16
#define cvac 2.99792458E8

/* --------- */
/* FUNCTIONS */
/* --------- */

void multiply(double complex mat1[Ndim][Ndim], double complex mat2[Ndim][Ndim], double complex res[Ndim][Ndim]);
void transpose(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim]);
void dagger(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim]);
void fancy(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim], int mode);
const double complex tr(double complex mat[Ndim][Ndim]);

double cquad(double my_f(), void * p, double a, double b, double relerr, double epsabs);
double qags(double my_f(), void * p, double a, double b, double relerr, double epsabs);
double integinf(double my_f(), void * p, double a, double relerr, double epsabs);
