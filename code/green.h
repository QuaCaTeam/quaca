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
const double complex hankalt0(double complex x);
const double complex hankalt1(double complex x);
const double complex hankaltn(int n, double complex x);
const double complex bessalt0(double complex x);
const double complex bessalt1(double complex x);
const double complex bessaltn(int n, double complex x);

const double complex epsilon(double complex omega, double omega_p, double gamma);
const double complex kappa(double complex omega, double c, double h);
const double complex eta(double complex omega, double omega_p, double gamma, double c, double h);

const double complex anum(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex bmm(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex bnn(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex bmn(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex cdenom(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex rMNn(int n, double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex refCoeffn(int sig, int n, double complex omega, double omega_p, double gamma, double c, double h, double R);

const double complex rNN0(double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex rNN1(double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex rMM1(double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex rMN1(double complex omega, double omega_p, double gamma, double c, double h, double R);
const double complex rNN0SF(double complex omega, double omega_p, double gamma, double c, double x);
const double complex rNN1SF(double complex omega, double omega_p, double gamma, double c, double x);

void greencent(double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double complex g[2][2]);
void greencentNF(double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double complex g[2][2]);

/* 
 * mathematical helping functions
 */
const double complex hankel1(double complex n, double complex x) {
    double complex result = ac_besselj(n,x) + I*ac_bessely(n,x);
    return result;
};

const double complex hankalt0(double complex x) {
    double complex result = -hankel1(1,x)/hankel1(0,x);
    return result;
};

const double complex hankalt1(double complex x) {
    double complex result = 1.0/2.0*(hankel1(0,x)-hankel1(2,x))/hankel1(1,x);
    return result;
};

const double complex hankaltn(int n, double complex x) {
    double complex result = (n/x*hankel1(n,x) - hankel1(n+1,x))/hankel1(n,x);
    return result;
};

const double complex bessalt0(double complex x) {
    double complex result = -ac_besselj(1,x)/ac_besselj(0,x);
    return result;
};

const double complex bessalt1(double complex x) {
    double complex result = 1.0/2.0*(ac_besselj(0,x)-ac_besselj(2,x))/ac_besselj(1,x);
    return result;
};

const double complex bessaltn(int n, double complex x) {
    double complex result = (n/x*ac_besselj(n,x) - ac_besselj(n+1,x))/ac_besselj(n,x);
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

const double complex rNN0(double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = -(hankel1(0,R*kappa(omega,c,h))*(hankalt0(R*kappa(omega,c,h))*bessalt0(R*kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h) - hankalt0(R*eta(omega,omega_p,gamma,c,h))*(hankalt0(R*kappa(omega,c,h)) + bessalt0(R*kappa(omega,c,h))*epsilon(omega,omega_p,gamma))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankalt0(R*eta(omega,omega_p,gamma,c,h))*hankalt0(R*eta(omega,omega_p,gamma,c,h))*epsilon(omega,omega_p,gamma)*kappa(omega,c,h)*kappa(omega,c,h)))/(ac_besselj(0,R*kappa(omega,c,h))*(bessalt0(R*kappa(omega,c,h))*bessalt0(R*kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h) - hankalt0(R*eta(omega,omega_p,gamma,c,h))*bessalt0(R*kappa(omega,c,h))*(1 + epsilon(omega,omega_p,gamma))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankalt0(R*eta(omega,omega_p,gamma,c,h))*hankalt0(R*eta(omega,omega_p,gamma,c,h))*epsilon(omega,omega_p,gamma)*kappa(omega,c,h)*kappa(omega,c,h)));
    return result;
};

const double complex rNN1(double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = -(hankel1(1,R*kappa(omega,c,h))*(hankalt1(R*kappa(omega,c,h))*bessalt1(R*kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h) - (h*h*omega*omega*(-1 + epsilon(omega,omega_p,gamma))*(-1 + epsilon(omega,omega_p,gamma)))/(c*c*R*R*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h)*kappa(omega,c,h)) - hankalt1(R*eta(omega,omega_p,gamma,c,h))*(hankalt1(R*kappa(omega,c,h)) + bessalt1(R*kappa(omega,c,h))*epsilon(omega,omega_p,gamma))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankalt1(R*eta(omega,omega_p,gamma,c,h))*hankalt1(R*eta(omega,omega_p,gamma,c,h))*epsilon(omega,omega_p,gamma)*kappa(omega,c,h)*kappa(omega,c,h)))/(ac_besselj(1,R*kappa(omega,c,h))*(bessalt1(R*kappa(omega,c,h))*bessalt1(R*kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h) - (h*h*omega*omega*(-1 + epsilon(omega,omega_p,gamma))*(-1 + epsilon(omega,omega_p,gamma)))/(c*c*R*R*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h)*kappa(omega,c,h)) - hankalt1(R*eta(omega,omega_p,gamma,c,h))*bessalt1(R*kappa(omega,c,h))*(1 + epsilon(omega,omega_p,gamma))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankalt1(R*eta(omega,omega_p,gamma,c,h))*hankalt1(R*eta(omega,omega_p,gamma,c,h))*epsilon(omega,omega_p,gamma)*kappa(omega,c,h)*kappa(omega,c,h)));
    return result;
};

const double complex rMM1(double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = -(hankel1(1,R*kappa(omega,c,h))*(hankalt1(R*kappa(omega,c,h))*bessalt1(R*kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h) - (h*h*omega*omega*(-1 + epsilon(omega,omega_p,gamma))*(-1 + epsilon(omega,omega_p,gamma)))/(c*c*R*R*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h)*kappa(omega,c,h)) - hankalt1(R*eta(omega,omega_p,gamma,c,h))*(bessalt1(R*kappa(omega,c,h)) + hankalt1(R*kappa(omega,c,h))*epsilon(omega,omega_p,gamma))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankalt1(R*eta(omega,omega_p,gamma,c,h))*hankalt1(R*eta(omega,omega_p,gamma,c,h))*epsilon(omega,omega_p,gamma)*kappa(omega,c,h)*kappa(omega,c,h)))/(ac_besselj(1,R*kappa(omega,c,h))*(bessalt1(R*kappa(omega,c,h))*bessalt1(R*kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h) - (h*h*omega*omega*(-1 + epsilon(omega,omega_p,gamma))*(-1 + epsilon(omega,omega_p,gamma)))/(c*c*R*R*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h)*kappa(omega,c,h)) - hankalt1(R*eta(omega,omega_p,gamma,c,h))*bessalt1(R*kappa(omega,c,h))*(1 + epsilon(omega,omega_p,gamma))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankalt1(R*eta(omega,omega_p,gamma,c,h))*hankalt1(R*eta(omega,omega_p,gamma,c,h))*epsilon(omega,omega_p,gamma)*kappa(omega,c,h)*kappa(omega,c,h)));
    return result;
};

const double complex rMN1(double complex omega, double omega_p, double gamma, double c, double h, double R) {
    double complex result = -(hankel1(1,R*kappa(omega,c,h))*I*h*omega/(R*c*kappa(omega,c,h))*(1-epsilon(omega,omega_p,gamma))*(hankalt1(R*eta(omega,omega_p,gamma,c,h))- bessalt1(R*eta(omega,omega_p,gamma,c,h))))/(ac_besselj(1,R*kappa(omega,c,h))*((bessalt1(R*kappa(omega,c,h))*bessalt1(R*kappa(omega,c,h))*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h) - (h*h*omega*omega*(-1 + epsilon(omega,omega_p,gamma))*(-1 + epsilon(omega,omega_p,gamma)))/(c*c*R*R*eta(omega,omega_p,gamma,c,h)*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h)*kappa(omega,c,h)) - hankalt1(R*eta(omega,omega_p,gamma,c,h))*bessalt1(R*kappa(omega,c,h))*(1 + epsilon(omega,omega_p,gamma))*eta(omega,omega_p,gamma,c,h)*kappa(omega,c,h) + hankalt1(R*eta(omega,omega_p,gamma,c,h))*hankalt1(R*eta(omega,omega_p,gamma,c,h))*epsilon(omega,omega_p,gamma)*kappa(omega,c,h)*kappa(omega,c,h))));
    return result;
};

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
void greencent(double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double complex g[2][2]) {

    // define non zero entries
    double complex gorth = I/(8*eps0)*(omega*omega/(c*c)*rMM1(omega,omega_p,gamma,c,h,R) + h*h*rNN1(omega,omega_p,gamma,c,h,R) - 2*I*omega*h/c*rMN1(omega,omega_p,gamma,c,h,R));
    double complex gzz = 2*kappa(omega,c,h)*kappa(omega,c,h)*rNN0(omega,omega_p,gamma,c,h,R);

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

// near field limit
void greencentNF(double complex omega, double omega_p, double gamma, double c, double h, double R, double eps0, double complex g[2][2]) {

    // define non zero entries
    double complex gorth = I*h*h/(8*eps0)*rNN1(omega,omega_p,gamma,c,h,R);
    double complex gzz = I*h*h/(8*eps0)*(-2)*rNN0(omega,omega_p,gamma,c,h,R);

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
