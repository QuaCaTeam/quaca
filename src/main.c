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
    double F0val, Fanarval, Fanatval, Ffreetval, Ffreerval;
    inputparams.eps0= 1./(4*PI);

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
    spac = pow(inputparams.stop/inputparams.start,1./inputparams.steps);

    /* Starting calculations */
    printf("\n===========================================\n");
    printf("CALCULATION STARTED!\n");

    clock_t c0 = clock();
    clock_t cl0, cl1;
    for (l=0; l<=inputparams.steps; ++l){

        inputparams.v = inputparams.start*pow(spac,l);

        /* Performing calculations */

        /* translational contribution */
        inputparams.transroll = 0;

        cl0 = clock();
        QFt = QF(IntQF, &inputparams);
        cl1 = clock();
        printf("Time elapsed QFt: %3.2f sec, ", (double)(cl1-cl0)/CLOCKS_PER_SEC);

        /* rolling contribution */
        inputparams.transroll = 1;

        cl0 = clock();
        QFr = QF(IntQF, &inputparams);
        cl1 = clock();
        printf("Time elapsed QFr: %3.2f sec\n", (double)(cl1-cl0)/CLOCKS_PER_SEC);

        /* analytics */ 
        F0val = F0(&inputparams);
        Fanatval = Fanat(&inputparams);
        Fanarval = Fanar(&inputparams);
        Ffreetval = Ffreet(&inputparams);
        Ffreerval = Ffreer(&inputparams);

        /* Print result to the screen */
        printf("v        | QFt/F0   | QFr/F0    | Fanat/F0 | Fanar/F0  | Ffreet/F0| Ffreer/F0\n");
        printf("%.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e\n\n",
                inputparams.v, QFt/F0val, QFr/F0val,Fanatval/F0val, Fanarval/F0val,
                Ffreetval/F0val, Ffreerval/F0val);

        /* write to the file */
        fprintf(fp, "%.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e\n",
                inputparams.v, QFt/F0val, QFr/F0val,Fanatval/F0val, Fanarval/F0val,
                Ffreetval/F0val, Ffreerval/F0val);
        fflush(fp);
    };

    clock_t c1 = clock();
    printf("\n------------------------------\n");

    /* Bye! */
    printf("\n");
    printf("Finished calculating %d points in %3.2f sec. \n", inputparams.steps, (c1-c0)/1.e6);
    printf("\nBYE!\n");
    printf("===========================================\n");

    return 0;
};

////////////////////////////////////////////////////////////////////////////////
