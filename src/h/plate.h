/* --------- */
/* LIBRARIES */
/* --------- */

#include "qfhelp.h"

/* universal constants */
double c, hbar, eps0, a0;

/* material and system parameters */
double wp1, wsp1, wa, gamMu, g1, einf, beta;
double vF, aF;
double v, za, QFt, QFr;

/* numeric parameters */
double relerr, recerr, abserr;
double kcut;

/* flags */
unsigned int transroll;
unsigned int muquest;
unsigned int omsgn;
unsigned int retard;

/* running variable */
char runvar[100];
char scale[100];
double start;
double stop;
unsigned int steps;


/* --------- */
/* FUNCTIONS */
/* --------- */

void input(char file[], unsigned int verbose);

void refl(double complex r[2], double w, double complex kap);
void reflhydro(double complex r[2] , double w, double complex kap);

double complex mu( double w);

void Gint(double complex Gten[Ndim][Ndim], double w, unsigned int RorI, unsigned int kx, unsigned int theta, unsigned int T);
void Gintnew(double complex Gten[Ndim][Ndim], double w, int RorI, int kx, int theta, int T);

void alpha(double complex alp[Ndim][Ndim], double w);

double IntQF(double w);
double IntQFnew( double w);
double IntQFfree( double w);

double QF(double IntQF());

/* analytics */
double F0(double wsp1, double a0, double eps0);
double Fanat(unsigned int muquest, double a0, double g1, double eps0, double wp1, double v, double za, double beta);
double Fanar(double a0, double g1, double eps0, double wp1, double v, double za, double beta);
double Ffreet(double a0, double g1, double eps0, double wp1, double v, double za, double beta);
double Ffreer(double a0, double g1, double eps0, double wp1, double v, double za, double beta);
