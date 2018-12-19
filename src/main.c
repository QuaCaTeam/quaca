////////////////////////////////////////////////////////////////////////////////
// This program is designed to calculate nonequilibrium quantum friction
// as well as the angular momentum and moment of inertia of a dipole hovering
// at a constartnt distartnce and velocity above a macroscopic surface described by
// a local dielectric function.
// For the calculation of the quantum friction we use the form:
//
// F_fric = \int_{0}^\infty dw (
//         -2Tr[S(w) \int d^2k/(2pi)^2 kx G_I(k,kx v + w)]
//         +2 hbar/pi Tr[alpI(w)\int d^2k/(2pi)^2 kx G_I(k,kx v + w)]]
//          *H(kx v + w) )
/******************************************************************************/
// The header.h contains all functions and subroutines needed for those
// calculations, as well as all needed parameters. Here a brief overview over
// the functions and subroutines contained:
//
// Finite integration      : integ(f,a,b,relerr)
// Matrix multiplication   : multiply(mat1,mat2,res)
// Trace of a matrix       : tr(mat)
// Complex-transposed      : dagger(mat,dagmat)
// 1/(2i)(A - A^+)         : fancyI(mat, matI)
// Reflection coefficient  : rR(w,k) and rI(w,k)
// Integrated Green tensor : Gint(Gten,w,RorI,T,kx)
// Polarizability          : alpha(alp,w)
#include "h/plate.h"

/******************************************************************************/
// Main program

int main (int argc, char *argv[]) {

    /* Plot parameters */
    double spac;                     // and the respective spacing

    /* Dummies */
    register unsigned int l;
    double QFt, QFr; // QFtfree, QFrfree;
    double F0, Fanar, Fanat, Ffreet, Ffreer;
    double eps0= 1./(4*PI);

    /* Greetings */
    printf("===========================================\n");
    printf("WELCOME TO QFNUM!\n");
    printf("===========================================\n");

    /* check input file */
    if (argc == 1) {
        printf("No file passed!\n");
        exit(0);
    } else {
        input(argv[1], 1);
    }

    /* create output file */
    char outfile[strlen(argv[1])+1];
    memset(outfile, '\0', sizeof(outfile)); //clear memory location
    strncpy(outfile, argv[1], strlen(argv[1])-3);
    strcat(outfile, ".out");

    FILE *fp;
    fp = fopen(outfile, "w");
    if (fp == NULL) {
        printf("Couldn't open output file \" %s \" for writing!\n", outfile);
        exit(0);
    };

    /* calculate spacing */
    spac = pow(stop/start,1./steps);

    /* Starting calculations */
    printf("\n===========================================\n");
    printf("CALCULATION STARTED!\n");
    printf("v               | QFt/F0          | QFr/F0           | Fanat/F0        | Fanar/F0         | Ffreet/F0       | Ffreer/F0\n");

    clock_t c0 = clock();
    for (l=0; l<=steps; ++l){

        v = start*pow(spac,l);

        /* Performing calculations */
        /* translational contribution */
        transroll = 0;

        clock_t cl0 = clock();
        QFt =  integ(IntQF,1E-5,0.9E0*wa,relerr, abserr);
        clock_t cl1 = clock();
        printf("Time elapsed for QFt (Part 1): %3.2f sec\n", (cl1-cl0)/1.e6);

        abserr = fabs(QFt)*1E-2;
        cl0 = clock();
        QFt += integ(IntQF,0.9E0*wa,wa,relerr, abserr);
        cl1 = clock();
        printf("Time elapsed for QFt (Part 2): %3.2f sec\n", (cl1-cl0)/1.e6);

        abserr = fabs(QFt)*1E-2;
        cl0 = clock();
        QFt += integ(IntQF,wa,wsp1,relerr, abserr);
        cl1 = clock();
        printf("Time elapsed for QFt (Part 3): %3.2f sec\n", (cl1-cl0)/1.e6);

        abserr = fabs(QFt)*1E-2;
        cl0 = clock();
        QFt += integinf(IntQF,wsp1,relerr, abserr);
        cl1 = clock();
        printf("Time elapsed for QFt (Part 4): %3.2f sec\n", (cl1-cl0)/1.e6);


        /* rolling contribution */
        abserr = 1E-200;
        transroll = 1;


        cl0 = clock();
        QFr =  integ(IntQF,0E0,0.9E0*wa,relerr, abserr);
        cl1 = clock();
        printf("Time elapsed for QFr (Part 1): %3.2f sec\n", (cl1-cl0)/1.e6);

        abserr = fabs(QFr)*1E-2;
        cl0 = clock();
        QFr += integ(IntQF,0.9E0*wa,wa,relerr, abserr);
        cl1 = clock();
        printf("Time elapsed for QFr (Part 2): %3.2f sec\n", (cl1-cl0)/1.e6);

        abserr = fabs(QFr)*1E-2;
        cl0 = clock();
        QFr += integ(IntQF,wa,wsp1,relerr, abserr);
        cl1 = clock();
        printf("Time elapsed for QFr (Part 3): %3.2f sec\n", (cl1-cl0)/1.e6);

        abserr = fabs(QFr)*1E-2;
        cl0 = clock();
        QFr += integinf(IntQF,wsp1,relerr, abserr);
        cl1 = clock();
        printf("Time elapsed for QFr (Part 4): %3.2f sec\n", (cl1-cl0)/1.e6);

        abserr = 1E-200;


        cl0 = clock();
        /* Calculate normalization constartnt */
        F0 = -3*pow(wsp1,5)*a0/(2*PI*eps0);

        /* Calculating analytical approximation for small velocities */
        Fanat = -63*a0*a0*pow(g1/(eps0*wp1*wp1),2)*pow(v,3)/(pow(PI,3)*pow(2*za,10));
        Fanar = -45*Fanat/63.;

        if (muquest == 1){
            Fanat += -45*creal(mu(0E0))*a0*g1*pow(v,3)*pow(2*za,-7)*4*PI/(pow(PI*wa*wp1,2));
            Fanat += -8*v*a0*g1*creal(mu(0E0))*4*PI/(pow(wp1*beta*wa,2)*pow(2*za,5));
        };

        Fanat = Fanat - a0*a0*pow(g1/(eps0*wp1*wp1),2)*6*v/(PI*beta*beta*pow(2*za,8));
        Fanar = Fanar + 3*a0*a0*pow(g1/(eps0*wp1*wp1),2)*v/(PI*beta*beta*pow(2*za,8));
        Ffreet= - a0*a0*pow(g1*4*PI/(eps0*wp1*wp1),2)*3*v/(PI*beta*beta*pow(2*za,8));
        Ffreer= -Ffreet/2.;

        cl1 = clock();
        printf("Time elapsed for F0, Fanat, Fanar, Ffreet, Ffreer: %3.5f sec\n", (cl1-cl0)/1.e6);

        /* Print result to the screen */
        printf("%.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e\n",
                v, QFt/F0, QFr/F0,Fanat/F0, Fanar/F0,
                Ffreet/F0, Ffreer/F0);
        /* write to the file */
        fprintf(fp, "%.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e\n",
                v, QFt/F0, QFr/F0,Fanat/F0, Fanar/F0,
                Ffreet/F0, Ffreer/F0);
        fflush(fp);
    };

    clock_t c1 = clock();
    printf("\n------------------------------\n");

    /* Bye! */
    printf("\n");
    printf("Finished calculating %d points in %3.2f sec. \n", steps, (c1-c0)/1.e6);
    printf("\nBYE!\n");
    printf("===========================================\n");

    return 0;
};

////////////////////////////////////////////////////////////////////////////////
