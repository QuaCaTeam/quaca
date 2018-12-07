/*!
 * \file qfhelp.h
 * \brief Common functions and parameters among the numerical calculation of quantum friction.
 * \author C. H. E. && M. O.
 *
 * Common functions.
 */

/* --------- */
/* Libraries */
/* --------- */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <string.h>


/* ---------- */
/* Parameters */
/* ---------- */

// mathematical constants

/*!
 * \def PI
 * \brief Constant Pi.
 */
#define PI 3.14159265358979323846

/*!
 * \def Ndim
 * \brief Matrix dimension.
 */
#define Ndim 3

/*!
 * \def PBSTR
 * \brief String for progress bar.
 *
 * \def PBWIDTH
 * \brief Width of progress bar.
 */
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


/* --------- */
/* Functions */
/* --------- */

/*!
 * \fn void printProg(double percentage);
 * \brief Prints progress bar.
 * \param percentage Percentage to display.
 * \return void
 */
void printProg(double percentage);

/*!
 * \fn void multiply(double complex mat1[Ndim][Ndim], double complex mat2[Ndim][Ndim], double complex res[Ndim][Ndim])
 * \brief Multiplies to square matrices of dimension Ndim.
 * \param mat1 First square matrix of dimension Ndim.
 * \param mat2 Second matrix of dimension Ndim.
 * \param res Resulting matrix of dimension Ndim.
 * \return void
 */
void multiply(double complex mat1[Ndim][Ndim], double complex mat2[Ndim][Ndim], double complex res[Ndim][Ndim]);

/*!
 * \fn void transpose(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim])
 * \brief Calculates the transpose of a square matrix of dimension Ndim.
 * \param mat Square matrix of dimension Ndim.
 * \param res Resulting matrix of dimension Ndim.
 * \return void
 */
void transpose(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim]);

/*!
 * \fn void dagger(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim])
 * \brief Calculates the conjugate transpose of a square matrix of dimension Ndim.
 * \param mat Square matrix of dimension Ndim.
 * \param res Resulting matrix of dimension Ndim.
 * \return void
 */
void dagger(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim]);

/*!
 * \fn void fancy(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim], int mode)
 * \brief Calculates fancy R or fancy I of a matrix.
 * \param mat Square matrix of dimension Ndim.
 * \param res Resulting matrix of dimension Ndim.
 * \param mode Option for R or I, mode = 1 equals R and mode = -1 equals I.
 * \return void
 */
void fancy(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim], int mode);

/*!
 * \fn const double complex tr(double complex mat[Ndim][Ndim])
 * \brief Computes the trace of a square matrix.
 * \param mat Square matrix of dimension Ndim.
 * \return Trace of the matrix.
 */
const double complex tr(double complex mat[Ndim][Ndim]);

/*!
 * \fn double integ(double my_f(), double a, double b, double relerr, double epsabs)
 * \brief Integration wrapper for finite integrations.
 * \param my_f() Function of one argument that shall be integrated.
 * \param a Lower integral bound.
 * \param b Upper integral bound.
 * \param relerr Relative error.
 * \param epsabs Absolute error.
 *
 * Uses the cquad routine of the gsl library. https://www.gnu.org/software/gsl/manual/html_node/CQUAD-doubly_002dadaptive-integration.html
 */
double integ(double my_f(), double a, double b, double relerr, double epsabs);

/*!
 * \fn double integinf(double my_f(), double a, double relerr, double epsabs)
 * \brief Integration wrapper for semi infinite integration.
 * \param my_f() Function of one argument that shall be integrated.
 * \param a Lower integral bound.
 * \param relerr Relative error.
 * \param epsabs Absolute error.
 */
double integinf(double my_f(), double a, double relerr, double epsabs);
