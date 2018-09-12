/*! 
 * \file cyl.h
 * \brief Header file for quantum friction in a cylinder.
 * \author C. H. Egerland
 * \bug No known bugs.
 *
 * Defines mathematical operations involving matrices, helping functions to define the reflection coefficients, reflection coefficients themselves and the Greens's tensor.
 */

/* --------- */
/* Libraries */
/* --------- */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf_bessel.h>
#include "arb.h"
#include "arbcmath.h"


/* ---------------- */
/* Fixed Parameters */
/* ---------------- */

// mathematical constants

/*! 
 * \def PI 
 * Constant Pi.
 */
#define PI 3.14159265358979323846

/*! 
 * \def Ndim 
 * Matrix dimension.
 */
#define Ndim 3

// material constants

/*!
 * \var omega_p
 * Plasma frequency in eV.
 *
 * \var gamma_p
 * Plasma resistivity in eV.
 */
double omega_p, gamma_p;

// universal constants

/*!
 * \var c
 * Speed of light in vacuum.
 *
 * \var eps0
 * Vacuum permittivity.
 */
double c, eps0;

// geometry

/*!
 * \var R
 * Radius of the cylinder.
 */
double R;


/* ------------------- */
/* Predefine functions */
/* ------------------- */

/*!
 * \fn void multiply(double complex mat1[Ndim][Ndim], double complex mat2[Ndim][Ndim], double complex res[Ndim][Ndim])
 * \brief Multiplies to square matrices of dimension Ndim.
 * \param mat1[Ndim][Ndim] First square matrix of dimension Ndim.
 * \param mat2[Ndim][Ndim] Second matrix of dimension Ndim.
 * \param res[Ndim][Ndim] Resulting matrix of dimension Ndim.
 * \return void
 */
void multiply(double complex mat1[Ndim][Ndim], double complex mat2[Ndim][Ndim], double complex res[Ndim][Ndim]);

/*!
 * \fn void dagger(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim])
 * \brief Calculates the conjugate transpose of a square matrix of dimension Ndim.
 * \param mat[Ndim][Ndim] Square matrix of dimension Ndim.
 * \param res[Ndim][Ndim] Resulting matrix of dimension Ndim.
 * \return void
 */
void dagger(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim]);

/*!
 * \fn void fancy(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim], int mode)
 * \brief Calculates fancy R or fancy I of a matrix.
 * \param mat[Ndim][Ndim] Square matrix of dimension Ndim.
 * \param res[Ndim][Ndim] Resulting matrix of dimension Ndim.
 * \param mode Option for R or I, mode = 1 equals R and mode = -1 equals I.
 * \return void
 */
void fancy(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim], int mode);

/*!
 * \fn const double complex tr(double complex mat[Ndim][Ndim])
 * \brief Computes the trace of a square matrix.
 * \param mat[Ndim][Ndim] Square matrix of dimension Ndim.
 * \return Trace of the matrix.
 */
const double complex tr(double complex mat[Ndim][Ndim]);

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
 * \fn const double complex rNN0SF(double complex omega, double h)
 * \brief Computes the reflection coefficient rNN0 in the small frequency limit.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex rNN0SF(double complex omega, double x);

/*!
 * \fn const double complex rNN1SF(double complex omega, double h)
 * \brief Computes the reflection coefficient rNN0 in the small frequency limit.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 */
const double complex rNN1SF(double complex omega, double x);

/*!
 * \fn void greenfull(int n, double complex omega, double h, double rho, double complex g[3][3])
 * \brief Computes the full Green's tensor of the cylinder.
 * \param n Order.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 * \param rho Coordinate rho (Radius).
 * \param g[Ndim][Ndim] Resulting tensor.
 */
void greenfull(int n, double complex omega, double h, double rho, double complex g[3][3]);

/*!
 * \fn void greencent(double complex omega, double h, double complex g[3][3])
 * \brief Computes the Green's tensor on the coaxial line of the cylinder.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 * \param g[Ndim][Ndim] Resulting tensor.
 */
void greencent(double complex omega, double h, double complex g[3][3]);

/*!
 * \fn void greencentNF(double complex omega, double h, double complex g[3][3])
 * \brief Computes the Green's tensor on the coaxial line of the cylinder in the near field limit.
 * \param omega Frequency.
 * \param h Wave vector (Fourier variable of z).
 * \param g[Ndim][Ndim] Resulting tensor.
 */
void greencentNF(double complex omega, double h, double complex g[3][3]);


/* --------- */
/* Functions */
/* --------- */

void multiply(double complex mat1[Ndim][Ndim], double complex mat2[Ndim][Ndim], double complex res[Ndim][Ndim]) {
    int i, j, k;
    for (i = 0; i < Ndim; i++) {
        for (j = 0; j < Ndim; j++) {
            res[i][j] = 0.;
            for (k = 0; k < Ndim; k++) {
                res[i][j] += mat1[i][k]*mat2[k][j]; 
            }
        }
    }

};

void dagger(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim]) {
    int i, j;
    for (i = 0; i < Ndim; i++) {
        for (j = 0; j < Ndim; j++) {
            res[j][i] = conj(mat[i][j]); 
        }
    }
};

void fancy(double complex mat[Ndim][Ndim], double complex res[Ndim][Ndim], int mode) {
    int i, j;
    for (i = 0; i < Ndim; i++) {
        for (j = 0; j < Ndim; j++) {
            res[i][j] = -0.5*I*( mat[i][j] + mode*conj(mat[j][i]) );
        }
    }
};

const double complex tr(double complex mat[Ndim][Ndim]) {
    int i;
    double complex res;
    res = 0.;
    for (i = 0; i < Ndim; i++) {
        res += mat[i][i]; 
    }
    return res;
};

const double complex hankel1(double complex n, double complex x) {
    double complex result = ac_besselj(n,x) + I*ac_bessely(n,x);
    return result;
};

const double complex hankaltn(int n, double complex x) {
    double complex result = (n/x*hankel1(n,x) - hankel1(n+1,x))/hankel1(n,x);
    return result;
};

const double complex bessaltn(int n, double complex x) {
    double complex result = (n/x*ac_besselj(n,x) - ac_besselj(n+1,x))/ac_besselj(n,x);
    return result;
};

const double complex ac_besselj_diff(int n, double complex x) {
    double complex result = n/x*ac_besselj(n,x) - ac_besselj(n+1,x);
    return result;
};

const double complex epsilon(double complex omega) {
    double complex result = 1.0 - omega_p*omega_p/(omega*omega + I*omega*gamma_p);
    return result;
};

const double complex kappa(double complex omega, double h) {
    double complex result = csqrt(omega*omega/(c*c) - h*h);
    return result;
};

const double complex eta(double complex omega, double h) {
    double complex result = csqrt(omega*omega/(c*c)*epsilon(omega) - h*h);
    return result;
};

const double complex anum(int n, double complex omega, double h) {
    double complex result = -n*n*omega*omega*h*h*pow(R,4)/(c*c)*(epsilon(omega)-1)*(epsilon(omega)-1);
    return result;
};

const double complex bmm(int n, double complex omega, double h) {
    double complex result = pow(R,6)*kappa(omega,h)*kappa(omega,h)*eta(omega,h)*eta(omega,h)*(epsilon(omega)*hankaltn(n,eta(omega,h))*hankaltn(n,eta(omega,h))*kappa(omega,h)*kappa(omega,h) - (hankaltn(n,eta(omega,h))*bessaltn(n,kappa(omega,h)) + hankaltn(n,eta(omega,h))*hankaltn(n,kappa(omega,h))*epsilon(omega))*eta(omega,h)*kappa(omega,h) + hankaltn(n,kappa(omega,h))*bessaltn(n,kappa(omega,h))*eta(omega,h)*eta(omega,h));
    return result;
};

const double complex bnn(int n, double complex omega, double h) {
    double complex result = pow(R,6)*kappa(omega,h)*kappa(omega,h)*eta(omega,h)*eta(omega,h)*(epsilon(omega)*hankaltn(n,eta(omega,h))*hankaltn(n,eta(omega,h))*kappa(omega,h)*kappa(omega,h) - (epsilon(omega)*hankaltn(n,eta(omega,h))*bessaltn(n,kappa(omega,h)) + hankaltn(n,eta(omega,h))*hankaltn(n,kappa(omega,h)))*eta(omega,h)*kappa(omega,h) + hankaltn(n,kappa(omega,h))*bessaltn(n,kappa(omega,h))*eta(omega,h)*eta(omega,h));
    return result;
};

const double complex bmn(int n, double complex omega, double h) {
    double complex result = -anum(n,omega,h) + I*h*n*omega*pow(R,5)/c*kappa(omega,h)*eta(omega,h)*eta(omega,h)*(1 - epsilon(omega))*(hankaltn(n,eta(omega,h)) - bessaltn(n,eta(omega,h))); 
    return result;
};

const double complex cdenom(int n, double complex omega, double h) {
    double complex result = pow(R,6)*kappa(omega,h)*kappa(omega,h)*eta(omega,h)*eta(omega,h)*(epsilon(omega)*hankaltn(n,eta(omega,h))*hankaltn(n,eta(omega,h))*kappa(omega,h)*kappa(omega,h) - (epsilon(omega)+1)*hankaltn(n,eta(omega,h))*bessaltn(n,kappa(omega,h))*eta(omega,h)*kappa(omega,h) + bessaltn(n,kappa(omega,h))*bessaltn(n,kappa(omega,h))*eta(omega,h)*eta(omega,h));
    return result;
};

const double complex refCoeffn(int sig, int n, double complex omega, double h) {
    double complex b; 

    switch (sig) {
        case 1:
            b = bmm(n,omega,h); 
            break;
        case 2:
            b = bnn(n,omega,h); 
            break;
        case 3:
            b = bmn(n,omega,h); 
            break;
        default:
            printf("Please choose a reflection coefficient!\n 1 = mm, 2 = nn, 3 = mn\n");
            break;
    };

    double complex result = - hankel1(n,kappa(omega,h)*R)/ac_besselj(n,kappa(omega,h)*R)*(anum(n,omega,h) + b)/(anum(n,omega,h) + cdenom(n,omega,h));
    return result; 
}; 

const double complex rNN0SF(double complex omega, double x) {
    double complex result = -(hankel1(0,I*x)/ac_besselj(0,I*x)) + (2*I*omega*gamma_p*hankel1(0,I*x))/(PI*omega_p*omega_p*x*ac_besselj(0,I*x)*ac_besselj(0,I*x)*hankel1(1,I*x));
    return result;
};

const double complex rNN1SF(double complex omega, double x) {
    double complex result = -(hankel1(1,I*x)/ac_besselj(1,I*x)) - (2*omega*gamma_p*hankel1(1,I*x))/(PI*omega_p*omega_p*ac_besselj(1,I*x)*ac_besselj(1,I*x)*(hankel1(1,I*x)-I*x*hankel1(0,I*x)));
    return result;
};

void greenfull(int N, double complex omega, double h, double rho, double complex g[3][3]) {

    // fill tensor with zeros
    int i,j;
    for (i = 0; i <= 2; i++) {
        for (j = 0; j <= 2; j++) {
            g[i][j] = 0;
        }
    }

    // define non zero entries
    double complex g11, g22, g33, g13, g31 = 0 + 0*I;

    // sum over n upto N
    int n;
    for (n = 0; n <= N; n++) {
        g11 += I/2*(n*n/(kappa(omega,h)*rho*kappa(omega,h)*rho)*ac_besselj(n,kappa(omega,h)*rho)*refCoeffn(1,n,omega,h) + c*c*h*h/(omega*omega)*ac_besselj_diff(n,kappa(omega,h)*rho)*ac_besselj_diff(n,kappa(omega,h)*rho)*refCoeffn(2,n,omega,h) - 2*I*c*h*n/(omega*rho*kappa(omega,h))*ac_besselj(n,kappa(omega,h)*rho)*ac_besselj_diff(n,kappa(omega,h)*rho)*refCoeffn(3,n,omega,h));

        g22 += I/2*(ac_besselj_diff(n,kappa(omega,h)*rho)*ac_besselj_diff(n,kappa(omega,h)*rho)*refCoeffn(1,n,omega,h) + c*c*h*h*n*n/(omega*omega*rho*rho*kappa(omega,h)*kappa(omega,h))*ac_besselj(n,kappa(omega,h)*rho)*ac_besselj(n,kappa(omega,h)*rho)*refCoeffn(2,n,omega,h) - 2*I*c*h*n/(omega*rho*kappa(omega,h))*ac_besselj(n,kappa(omega,h)*rho)*ac_besselj_diff(n,kappa(omega,h)*rho)*refCoeffn(3,n,omega,h));
        
        g33 += I/2*(c*c/(omega*omega)*kappa(omega,h)*kappa(omega,h)*ac_besselj(n,kappa(omega,h)*rho)*ac_besselj(n,kappa(omega,h)*rho)*refCoeffn(2,n,omega,h));
        
        // prime in the sum means the 0 term has to be multiplied by 1/2 
        if (n == 0) {
            g11 = 0.5*g11; 
            g22 = 0.5*g22; 
            g33 = 0.5*g33; 
            g13 = 0.5*g13; 
            g31 = 0.5*g31; 
        }
    }

    // fill in non zero entries and add additional factor
    g[0][0] = omega*omega/(c*c*eps0)*g11;
    g[1][1] = omega*omega/(c*c*eps0)*g22;
    g[2][2] = omega*omega/(c*c*eps0)*g33;
    g[0][2] = omega*omega/(c*c*eps0)*g13;
    g[2][0] = omega*omega/(c*c*eps0)*g31;
};

void greencent(double complex omega, double h, double complex g[3][3]) {

    // define non zero entries
    double complex gorth = I/(8*eps0)*(omega*omega/(c*c)*refCoeffn(1,1,omega,h) + h*h*refCoeffn(2,1,omega,h) - 2*I*omega*h/c*refCoeffn(3,1,omega,h));
    double complex gzz = I/(8*eps0)*2*kappa(omega,h)*kappa(omega,h)*refCoeffn(2,0,omega,h);

    // fill tensor with zeros
    int i,j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            g[i][j] = 0;
        }
    };

    // fill in non zero entries
    g[0][0] = gorth;
    g[1][1] = gorth;
    g[2][2] = gzz;
};

// near field limit YOU FORGOT TO NEAR FIELD THE REFLECTION COEFFICIENTS
void greencentNF(double complex omega, double h, double complex g[3][3]) {

    // define non zero entries
    double complex gorth = I*h*h/(8*eps0)*refCoeffn(2,1,omega,h);
    double complex gzz = I*h*h/(8*eps0)*(-2)*refCoeffn(2,0,omega,h);

    // fill tensor with zeros
    int i,j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            g[i][j] = 0;
        }
    };

    // fill in non zero entries
    g[0][0] = gorth;
    g[1][1] = gorth;
    g[2][2] = gzz;
};
