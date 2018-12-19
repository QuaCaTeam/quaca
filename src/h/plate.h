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

void input(char file[], int verbose);

void refl(double complex r[2], double w, double complex kap);
void reflhydro(double complex r[2] , double w, double complex kap);

double complex mu( double w);

void Gint(double complex Gten[Ndim][Ndim], double w, int RorI, int kx, int theta, int T);
void Gintnew(double complex Gten[Ndim][Ndim], double w, int RorI, int kx, int theta, int T);

void alpha(double complex alp[Ndim][Ndim], double w);

double IntQF(double w);
double IntQFnew( double w);
double IntQFfree( double w);
