/*! 
 * \file cyl.c
 * \brief Functions for quantum friction in a cylinder.
 * \author C. H. Egerland
 * \bug greencentNF not right, reflection coefficients used are NOT in near field limit.
 *
 * Functions, for more explicit documentation see cyl.h
 */
 
/* --------- */
/* Libraries */
/* --------- */

#include "cyl.h"
#include "qfhelp.h"


/* --------- */
/* Functions */
/* --------- */

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
            exit(0);
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

void greenfull(int N, double complex omega, double h, double rho, double complex g[Ndim][Ndim]) {

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

void greencent(double complex omega, double h, double complex g[Ndim][Ndim]) {

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

void greencentNF(double complex omega, double h, double complex g[Ndim][Ndim]) {

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

void Greenint(double complex gres[Ndim][Ndim], double omega, int RorI, int horNoh, int theta) {
    
    // dummies
    double lim1, lim2;
    int i, j, sig;
    double complex gint[Ndim][Ndim], g[Ndim][Ndim];
    
    // integration limits
    lim1 = 0.;
    lim2 = hcut/(2*R);

    double wraph (double h) {
    
        // dummy integrand tensor
        double resh;
        int i, j;

        // initialize Green's tensor in g
        greencent(omega, h, g);

        /* Modify g -> gint */
        // make Green's tensor fancy R or I
        if (RorI == 0) {
            fancy(g, gint, 1);
        } else if (RorI == 1) {
            fancy(g, gint, -1); 
        } else {
            printf("Wrong RorI!\n");
            exit(0);
        }
        
        // add additional prefactor h in front
        if (horNoh == 0) {
        } else if (horNoh == 1) {
            for (i = 0; i < Ndim; i++) {
                for (j = 0; j < Ndim; j++) {
                    gint[i][j] = h*gint[i][j]; 
                }
            }
        } else {
            printf("Wrong horNoh!\n");
            exit(0);
        } 

        // add theta function
        if (theta == 0) {
        } else if (theta == 1) {
        } else {
            printf("Wrong theta!\n");
            exit(0);
        }

        // return the indice of interest
        if (sig == 1) {
            resh = gint[0][0];
        } else if (sig == 2) {
            resh = gint[1][1];
        } else if (sig == 3) {
            resh = gint[2][2];
        } else {
            printf("Wrong sig!\n");
            exit(0);
        } 

        return resh;
    };

    // flush tensor
    for (i = 0; i < Ndim; i++) {
        for (j = 0; j < Ndim; j++) {
            gres[i][j] = 0;
        }
    };

    // non zero entries
    sig = 1;
    gres[0][0] = integ(wraph, lim1, lim2, relerr);
    sig = 2;
    gres[1][1] = integ(wraph, lim1, lim2, relerr);
    sig = 3;
    gres[2][2] = integ(wraph, lim1, lim2, relerr);
};
