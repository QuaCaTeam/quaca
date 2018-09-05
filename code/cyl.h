// include libraries
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf_bessel.h>
#include "arb.h"
#include "arbcmath.h"

// constants
#define PI 3.14159265358979323846

// predefine functions
const double complex hankel1(double complex n, double complex x);
const double complex hankaltn(int n, double complex x);
const double complex bessaltn(int n, double complex x);
const double complex ac_besselj_diff(int n, double complex x);

const double complex epsilon(double complex omega, double omega_p, double gamma);
const double complex kappa(double complex omega, double c, double h);
const double complex eta(double complex omega, double omega_p, double gamma, double c, double h);

const double complex anum(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex bmm(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex bnn(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex bmn(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex cdenom(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex refCoeffn(int sig, int n, double complex omega, double omega_p, double gamma, double c, double h, double R);

const double complex rNN0SF(double complex omega, double omega_p, double gamma, double c, double x);
const double complex rNN1SF(double complex omega, double omega_p, double gamma, double c, double x);

void greenfull(int N, double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double rho, double complex g[2][2]);
void greencent(double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double complex g[2][2]);
void greencentNF(double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double complex g[2][2]);

/* 
 * mathematical helping functions
 */
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

/*
 * physical helping functions
 */

const double complex epsilon(double complex omega, double omega_p, double gamma) {
    double complex result = 1.0 - omega_p*omega_p/(omega*omega + I*omega*gamma);
    return result;
};

const double complex kappa(double complex omega, double c, double h) {
    double complex result = csqrt(omega*omega/(c*c) - h*h);
    return result;
};

const double complex eta(double complex omega, double omega_p, double gamma, double c, double h) {
    double complex result = csqrt(omega*omega/(c*c)*epsilon(omega,omega_p,gamma) - h*h);
    return result;
};


/*
 * general reflection coefficients
 */

const double complex anum(int n, double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = -n*n*omega*omega*h*h*pow(R,4)/(c*c)*(epsilon(omega, omega_p, gamma)-1)*(epsilon(omega, omega_p, gamma)-1);
    return result;
};

const double complex bmm(int n, double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = pow(R,6)*kappa(omega,c,h)*kappa(omega,c,h)*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*(epsilon(omega,omega_p,gamma)*hankaltn(n,eta(omega,omega_p,gamma,c,h))*hankaltn(n,eta(omega,omega_p,gamma,c,h))*kappa(omega,c,h)*kappa(omega,c,h) - (hankaltn(n,eta(omega,omega_p,gamma,c,h))*bessaltn(n,kappa(omega,c,h)) + hankaltn(n,eta(omega,omega_p,gamma,c,h))*hankaltn(n,kappa(omega,c,h))*epsilon(omega,omega_p,gamma))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankaltn(n,kappa(omega,c,h))*bessaltn(n,kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h));
    return result;
};

const double complex bnn(int n, double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = pow(R,6)*kappa(omega,c,h)*kappa(omega,c,h)*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*(epsilon(omega,omega_p,gamma)*hankaltn(n,eta(omega,omega_p,gamma,c,h))*hankaltn(n,eta(omega,omega_p,gamma,c,h))*kappa(omega,c,h)*kappa(omega,c,h) - (epsilon(omega,omega_p,gamma)*hankaltn(n,eta(omega,omega_p,gamma,c,h))*bessaltn(n,kappa(omega,c,h)) + hankaltn(n,eta(omega,omega_p,gamma,c,h))*hankaltn(n,kappa(omega,c,h)))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankaltn(n,kappa(omega,c,h))*bessaltn(n,kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h));
    return result;
};

const double complex bmn(int n, double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = -anum(n,omega,omega_p,gamma,c,h,R) + I*h*n*omega*pow(R,5)/c*kappa(omega,c,h)*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*(1 - epsilon(omega,omega_p,gamma))*(hankaltn(n,eta(omega,omega_p,gamma,c,h)) - bessaltn(n,eta(omega,omega_p,gamma,c,h))); 
    return result;
};

const double complex cdenom(int n, double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = pow(R,6)*kappa(omega,c,h)*kappa(omega,c,h)*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*(epsilon(omega,omega_p,gamma)*hankaltn(n,eta(omega,omega_p,gamma,c,h))*hankaltn(n,eta(omega,omega_p,gamma,c,h))*kappa(omega,c,h)*kappa(omega,c,h) - (epsilon(omega,omega_p,gamma)+1)*hankaltn(n,eta(omega,omega_p,gamma,c,h))*bessaltn(n,kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + bessaltn(n,kappa(omega,c,h))*bessaltn(n,kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h));
    return result;
};

const double complex refCoeffn(int sig, int n, double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex b; 

    switch (sig) {
        case 1:
            b = bmm(n,omega,omega_p,gamma,c,h,R); 
            break;
        case 2:
            b = bnn(n,omega,omega_p,gamma,c,h,R); 
            break;
        case 3:
            b = bmn(n,omega,omega_p,gamma,c,h,R); 
            break;
        default:
            printf("Please choose a reflection coefficient!\n 1 = mm, 2 = nn, 3 = mn\n");
            break;
    };

    double complex result = - hankel1(n,kappa(omega,c,h)*R)/ac_besselj(n,kappa(omega,c,h)*R)*(anum(n,omega,omega_p,gamma,c,h,R) + b)/(anum(n,omega,omega_p,gamma,c,h,R) + cdenom(n,omega,omega_p,gamma,c,h,R));
    return result; 
}; 


/*
 * specific reflection coefficients
 */

const double complex rNN0SF(double complex omega, double omega_p, double gamma, double c, double x) {
    double complex result = -(hankel1(0,I*x)/ac_besselj(0,I*x)) + (2*I*omega*gamma*hankel1(0,I*x))/(PI*omega_p*omega_p*x*ac_besselj(0,I*x)*ac_besselj(0,I*x)*hankel1(1,I*x));
    return result;
};

const double complex rNN1SF(double complex omega, double omega_p, double gamma, double c, double x) {
    double complex result = -(hankel1(1,I*x)/ac_besselj(1,I*x)) - (2*omega*gamma*hankel1(1,I*x))/(PI*omega_p*omega_p*ac_besselj(1,I*x)*ac_besselj(1,I*x)*(hankel1(1,I*x)-I*x*hankel1(0,I*x)));
    return result;
};


/*
 * Green's tensor
 */

void greenfull(int N, double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double rho, double complex g[2][2]) {

    // fill tensor with zeros
    int i,j;
    for (i = 0; i <= 2; i++) {
        for (j = 0; j <= 2; j++) {
            g[i][j] = 0;
        }
    };

    // define non zero entries
    double complex g11, g22, g33, g13, g31 = 0 + 0*I;

    // sum over n upto N
    int n;
    for (n = 0; n <= N; n++) {
        g11 += I/2*(n*n/(kappa(omega,c,h)*rho*kappa(omega,c,h)*rho)*ac_besselj(n,kappa(omega,c,h)*rho)*refCoeffn(1,n,omega,omega_p,gamma,c,h,R) + c*c*h*h/(omega*omega)*ac_besselj_diff(n,kappa(omega,c,h)*rho)*ac_besselj_diff(n,kappa(omega,c,h)*rho)*refCoeffn(2,n,omega,omega_p,gamma,c,h,R) - 2*I*c*h*n/(omega*rho*kappa(omega,c,h))*ac_besselj(n,kappa(omega,c,h)*rho)*ac_besselj_diff(n,kappa(omega,c,h)*rho)*refCoeffn(3,n,omega,omega_p,gamma,c,h,R));

        g22 += I/2*(ac_besselj_diff(n,kappa(omega,c,h)*rho)*ac_besselj_diff(n,kappa(omega,c,h)*rho)*refCoeffn(1,n,omega,omega_p,gamma,c,h,R) + c*c*h*h*n*n/(omega*omega*rho*rho*kappa(omega,c,h)*kappa(omega,c,h))*ac_besselj(n,kappa(omega,c,h)*rho)*ac_besselj(n,kappa(omega,c,h)*rho)*refCoeffn(2,n,omega,omega_p,gamma,c,h,R) - 2*I*c*h*n/(omega*rho*kappa(omega,c,h))*ac_besselj(n,kappa(omega,c,h)*rho)*ac_besselj_diff(n,kappa(omega,c,h)*rho)*refCoeffn(3,n,omega,omega_p,gamma,c,h,R));
        
        g33 += I/2*(c*c/(omega*omega)*kappa(omega,c,h)*kappa(omega,c,h)*ac_besselj(n,kappa(omega,c,h)*rho)*ac_besselj(n,kappa(omega,c,h)*rho)*refCoeffn(2,n,omega,omega_p,gamma,c,h,R));
        
        // prime in the sum means the 0 term has to be multiplied by 1/2 
        if (n == 0) {
            g11 = 0.5*g11; 
            g22 = 0.5*g22; 
            g33 = 0.5*g33; 
            g13 = 0.5*g13; 
            g31 = 0.5*g31; 
        }
    };

    // fill in non zero entries and add additional factor
    g[0][0] = omega*omega/(c*c*eps0)*g11;
    g[1][1] = omega*omega/(c*c*eps0)*g22;
    g[2][2] = omega*omega/(c*c*eps0)*g33;
    g[0][2] = omega*omega/(c*c*eps0)*g13;
    g[2][0] = omega*omega/(c*c*eps0)*g31;
};

void greencent(double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double complex g[2][2]) {

    // define non zero entries
    double complex gorth = I/(8*eps0)*(omega*omega/(c*c)*refCoeffn(1,1,omega,omega_p,gamma,c,h,R) + h*h*refCoeffn(2,1,omega,omega_p,gamma,c,h,R) - 2*I*omega*h/c*refCoeffn(3,1,omega,omega_p,gamma,c,h,R));
    double complex gzz = I/(8*eps0)*2*kappa(omega,c,h)*kappa(omega,c,h)*refCoeffn(2,0,omega,omega_p,gamma,c,h,R);

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
void greencentNF(double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double complex g[2][2]) {

    // define non zero entries
    double complex gorth = I*h*h/(8*eps0)*refCoeffn(2,1,omega,omega_p,gamma,c,h,R);
    double complex gzz = I*h*h/(8*eps0)*(-2)*refCoeffn(2,0,omega,omega_p,gamma,c,h,R);

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
