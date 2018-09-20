/*! 
 * \file cyl.h
 * \brief Header file for quantum friction in a cylinder.
 * \author C. H. Egerland
 * \bug greencentNF is still not right, since the reflection coefficients used are not in the NEAR FIELD limit.
 *
 * Defines mathematical operations involving matrices, helping functions to define the reflection coefficients, reflection coefficients themselves and the Greens's tensor.
 */

/* --------- */
/* Libraries */
/* --------- */

#include "qfhelp.h"

#include <gsl/gsl_sf_bessel.h>
#include "arb.h"
#include "arbcmath.h"


// material constants

/*!
 * \var omega_p
 * \brief Plasma frequency in eV.
 *
 * \var gamma_p
 * \brief Plasma resistivity in eV.
 */
double omega_p, gamma_p;

// universal constants

/*!
 * \var c
 * \brief Speed of light in vacuum.
 *
 * \var eps0
 * \brief Vacuum permittivity.
 *
 * \var hbar
 * \brief Reduced Plancks constant.
 */
double c, eps0, hbar;

// miscellaneous

/*!
 * \var R
 * \brief Radius of the cylinder.
 */
double R;

/*!
 * \var hcut
 * \brief Cutoff for h integration.
 */
double hcut;

/*!
 * \var relerr
 * \brief Relative error of the integration routine.
 */
double relerr, absr;

/* ------------------- */
/* Predefine functions */
/* ------------------- */

/*!
 * \fn const double complex hankel1(double complex n, double complex x)
 * \brief Computes the Hankel function of first kind.
 * \param n Order.
 * \param x Argument.
 */
const double complex hankel1(double complex n, double complex x);

/*!
 * \fn const double complex hankaltn(int n, double complex x)
 * \brief Computes derivative of hankel1 devided by hankel1.
 * \param n Order.
 * \param x Argument.
 */
const double complex hankaltn(int n, double complex x);

/*!
 * \fn const double complex bessaltn(int n, double complex x)
 * \brief Computes derivative of besselJ devided by besselJ.
 * \param n Order.
 * \param x Argument.
 */
const double complex bessaltn(int n, double complex x);

/*!
 * \fn const double complex ac_besselj_diff(int n, double complex x)
 * \brief Computes derivative of besselJ.
 * \param n Order.
 * \param x Argument.
 */
const double complex ac_besselj_diff(int n, double complex x);

/*!
 * \fn const double complex epsilon(double complex omega)
 * \brief Computes the permittivity resulting from the Drude model.
 * \param omega Frequency.
 */
const double complex epsilon(double complex omega);

/*!
 * \fn const double complex kappa(double complex omega, double h)
 * \brief Computes kappa.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex kappa(double complex omega, double h);

/*!
 * \fn const double complex eta(double complex omega, double h)
 * \brief Computes eta.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex eta(double complex omega, double h);

/*!
 * \fn const double complex anum(int n, double complex omega, double h)
 * \brief Computes helping function for reflection coefficients.
 * \param n Order.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex anum(int n, double complex omega, double h);

/*!
 * \fn const double complex bmm(int n, double complex omega, double h)
 * \brief Computes helping function for reflection coefficients.
 * \param n Order.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex bmm(int n, double complex omega, double h);

/*!
 * \fn const double complex bnn(int n, double complex omega, double h)
 * \brief Computes helping function for reflection coefficients.
 * \param n Order.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex bnn(int n, double complex omega, double h);

/*!
 * \fn const double complex bmn(int n, double complex omega, double h)
 * \brief Computes helping function for reflection coefficients.
 * \param n Order.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex bmn(int n, double complex omega, double h);

/*!
 * \fn const double complex cdenom(int n, double complex omega, double h)
 * \brief Computes helping function for reflection coefficients.
 * \param n Order.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex cdenom(int n, double complex omega, double h);

/*!
 * \fn const double complex refCoeffn(int sig, int n, double complex omega, double h)
 * \brief Computes the reflection coefficients in a cylinder.
 * \param sig Mode.
 * \param n Order.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 *
 * Choose sig = 1 for MM, sig = 2 for NN or sig = 3 for MN.
 */
const double complex refCoeffn(int sig, int n, double complex omega, double h);

/*!
 * \fn const double complex rNN0SF(double complex omega, double x)
 * \brief Computes the reflection coefficient rNN0 in the small frequency limit.
 * \param omega Frequency.
 * \param x Dimensionless wave vector (x = hR).
 */
const double complex rNN0SF(double complex omega, double x);

/*!
 * \fn const double complex rNN1SF(double complex omega, double x)
 * \brief Computes the reflection coefficient rNN0 in the small frequency limit.
 * \param omega Frequency.
 * \param x Dimensionless wave vector (x = hR).
 */
const double complex rNN1SF(double complex omega, double x);

/*!
 * \fn void greenfull(int N, double complex omega, double h, double rho, double complex g[Ndim][Ndim])
 * \brief Computes the full Green's tensor of the cylinder.
 * \param N Order.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 * \param rho Coordinate rho (Radius).
 * \param g Resulting tensor.
 */
void greenfull(int N, double complex omega, double h, double rho, double complex g[Ndim][Ndim]);

/*!
 * \fn void greencent(double complex omega, double h, double complex g[Ndim][Ndim])
 * \brief Computes the Green's tensor on the coaxial line of the cylinder.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 * \param g Resulting tensor.
 */
void greencent(double complex omega, double h, double complex g[Ndim][Ndim]);

/*!
 * \fn void greencentNF(double complex omega, double h, double complex g[Ndim][Ndim])
 * \brief Computes the Green's tensor on the coaxial line of the cylinder in the near field limit.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 * \param g Resulting tensor.
 */
void greencentNF(double complex omega, double h, double complex g[Ndim][Ndim]);

/*!
 * \fn void Greenint(double complex gres[Ndim][Ndim], double omega, int RorI, int horNoh, int theta)
 * \brief Computes the integral of greencent over h.
 * \param gres Resulting tensor.
 * \param omega Frequency.
 * \param RorI Flag for fancy R or I.
 * \param horNoh Flag for additional prefactor h.
 * \param theta Flag for theta function.
 */
void Greenint(double complex gres[Ndim][Ndim], double omega, int RorI, int horNoh, int theta);
