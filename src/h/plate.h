/* --------- */
/* LIBRARIES */
/* --------- */

#include "qfhelp.h"

struct parameters {
    /* universal constants */
    double eps0, a0;

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
};

struct parameters inputparams;

/* --------- */
/* FUNCTIONS */
/* --------- */

void input(char file[], unsigned int verbose);

void refl(double complex r[2], double w, double complex kap, void * p);

void Gint(double complex Gten[Ndim][Ndim], double w, void * p, unsigned int RorI, unsigned int kx, unsigned int theta, unsigned int T);

void alpha(double complex alp[Ndim][Ndim], double w, void * p);
double complex mu(double w, void * p);

double IntQF(double w, void * p);

double QF(double IntQF(), void * p);

/* analytics */
double F0(void * p);
double Fanat(void * p);
double Fanar(void * p);
double Ffreet(void * p);
double Ffreer(void * p);
