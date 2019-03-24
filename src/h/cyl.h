/* --------- */
/* LIBRARIES */
/* --------- */

#include "qfhelp.h"
#include "arbcmath.h"

struct parametersCyl {
    /* universal constants */
    double eps0, a0;

    /* material and system parameters */
    double wp1, wsp1, wa, gamMu, g1, einf, beta;
    double vF, aF;
    double v, R, QFt, QFr;

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

struct parametersCyl inputparamsCyl;


/* --------- */
/* FUNCTIONS */
/* --------- */

void inputCyl(char file[], unsigned int verbose);
const double complex hankel1(double complex n, double complex x);
const double complex hn(int n, double complex x);
const double complex besseljn(int n, double complex x);


double complex rNN(int n, double w, double h, void * p);
double complex rMM(int n, double w, double h, void * p);
double complex rMN(int n, double w, double h, void * p);

double complex rNNNF(int n, double w, double h, void * p);

void GCNFint(double Gten[Ndim][Ndim], double w, void * p, int RorI);
void GCint(double Gten[Ndim][Ndim], double w, void * p, int RorI);

