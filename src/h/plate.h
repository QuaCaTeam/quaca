/*!
 * \file plate.h
 * \brief Header file for quantum friction of a plate.
 * \author M. O.
 */

/* --------- */
/* Libraries */
/* --------- */

#include "qfhelp.h"

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
double c, hbar, eps0, a0;

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
 * \var bet1
 * \brief sound velocity of the hydrodynamic model proportional to vF.
 *
 * \var beta
 * \brief Temperature in eV.
 *
 * \var einf
 * \brief Background permittivity.
 */
double wp1, wsp1, wa, g1, bet1, einf, beta;

/*!
 * \var vF
 * \brief Fermi velocity.
 *
 * \var aF
 * \brief Fine structure constant.
 *
 * \var Delt
 * \brief Band gap.
 *
 * \var OmD
 * \brief Respective band gap frequency.
 */
double vF, aF, Delt, OmD;

/*!
 * \var kcut
 * \brief Cut-off for k integration.
 */
double  kcut;

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
double v, za, QFt, QFr;

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
double relerr, recerr, abserr;

/*!
 * \var transroll
 * \brief Flag for translational or rolling part.
 */
int transroll;


/* --------- */
/* Functions */
/* --------- */

/*!
 * \fn void input(char file[]);
 * \brief Parses file for variables.
 * \param file File to parse.
 * \return void
 */
void input();

/*!
 * \fn void refl(double complex r[2] ,double k, double w, double complex kap)
 * \brief
 * \param r
 * \param w
 * \param kap
 * \return
 */
void refl(double complex r[2], double w, double complex kap);

/*!
 * \fn void refl(double complex r[2] ,double k, double w, double complex kap)
 * \brief
 * \param r
 * \param w
 * \param kap
 * \return
 */
void reflhydro(double complex r[2] , double w, double complex kap);

/*!
 * \fn void reflGraph(double complex r[2] ,double k, double w, double complex kap)
 * \brief
 * \param r
 * \param w
 * \param kap
 * \return
 */
void reflGraph(double complex r[2] , double w, double complex kap);

/*!
 * \fn void Gint(double complex Gten[Ndim][Ndim], double w, int RorI, int kx, int theta, int T)
 * \brief
 * \param Gten
 * \param w
 * \param RorI
 * \param kx
 * \param theta
 * \param T
 * \return
 */
void Gint(double complex Gten[Ndim][Ndim], double w, int RorI, int kx, int theta, int T);

/*!
 * \fn void alpha(double complex alp[Ndim][Ndim], double w)
 * \brief
 * \param alp
 * \param w
 * \return
 */
void alpha(double complex alp[Ndim][Ndim], double w);

/*!
 * \fn double IntQF(double w)
 * \brief
 * \param w
 * \return
 */
double IntQF(double w);

/*!
 * \fn double IntQFfree(double w)
 * \brief
 * \param w
 * \return
 */
double IntQFfree( double w);