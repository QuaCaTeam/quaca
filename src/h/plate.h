/*! 
 * \file plate.h
 * \brief Header file for quantum friction of a plate.
 * \author M. O.
 *
 */

/* --------- */
/* Libraries */
/* --------- */

#include "qfhelp.h"


/* ---------- */
/* Parameters */
/* ---------- */


extern double complex Ly[3][3] = {
   {0.   , 0., 1.*_Complex_I} ,   /*  initializers for row indexed by 0 */
   {0.   , 0., 0.  } ,   /*  initializers for row indexed by 1 */
   {-1.*_Complex_I, 0., 0.  }     /*  initializers for row indexed by 2 */
};
extern double complex Ly2[3][3] = {
   {1.   , 0., 0.  } ,   /*  initializers for row indexed by 0 */
   {0.   , 0., 0.  } ,   /*  initializers for row indexed by 1 */
   {0.   , 0., 1.  }     /*  initializers for row indexed by 2 */
};


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
 *
 * \var a0
 * \brief Static polarizability.
 */
extern double c, hbar, eps0, a0;

// material and system parameters

/*!
 * \var wp1 
 * \brief Plasma frequency in eV.
 *
 * \var wsp1 
 * \brief Surface plasmon frequency in eV.
 *
 * \var wa
 * \brief Dipole resonance frequency in eV.
 *
 * \var g1
 * \brief Damping of material in eV.
 *
 * \var beta
 * \brief Temperature in eV.
 *
 * \var delta
 * \brief I don't know.
 *
 * \var einf
 * \brief Background permittivity.
 */
extern double wp1, wsp1, wa, g1, beta, delta, einf;

/*!
 * \var kcut
 * \brief Cut-off for k integration.
 */
extern double  kcut;

/*!
 * \var v
 * \brief Velocity in units of c.
 *
 * \var za
 * \brief Height of dipole in 1/eV
 *
 * \var QFt
 * \brief dont know.
 *
 * \var QFr
 * \brief dont know.
 */
extern double v, za, QFt, QFr;

/*!
 * \var relerr
 * \brief Relative error.
 *
 * \var recerr
 * \brief Increase of relerr per layer of integration.
 *
 * \var absr
 * \brief Absolute error.
 */
extern double relerr, recerr, absr=1e-200;

/*! 
 * \var transroll
 * \brief Flag for translational or rolling part.
 */
int transroll;


/* --------- */
/* Functions */
/* --------- */

double rI (double w, double k);
double rR (double w, double k);

void Gint(double complex Gten[Ndim][Ndim], double w, int RorI, int kx, int theta, int T);
void alpha(double complex alp[Ndim][Ndim], double w);
double anaAngL(double v);
double AngL(double w);
double Iner(double w);
double IntQF(double w);
