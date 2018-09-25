////////////////////////////////////////////////////////////////////////////////
// This program is designed to calculate nonequilibrium quantum friction
// as well as the angular momentum and moment of inertia of a dipole hovering
// at a constant distance and velocity above a macroscopic surface described by
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
int main () {

/* Plot parameters */
int maxi=1000;                        // plot points
double sta = 1.29999953 ;                 // start value of the calculation,
double sto = 1.299999533;                 // final value of the calculation
double spac =(sto-sta)/maxi;      // and the respective spacing

/* Dummies */
int l;
double w;
double prog;

/* Greetings */
printf("===========================================\n");
printf("QFNUM STARTED!\n\n");

/* System parameters (input routine is not implemented yet) */
input("../data/input/paratest.dat");

/* open the file */
FILE *fp; // output file
fp = fopen("../data/output/paratestout.dat", "w");
if (fp == NULL) {
    printf("I couldn't open results.dat for writing.\n");
    exit(0);
}


/* Starting calculations */
printf("\n------------------------------\n");
printf("STARTING CALCULATION\n");

clock_t c0 = clock();
for (l=0; l<=maxi; ++l){
    /* Evaluation */
    w = sta+spac*l;
    fprintf(fp, "%.10e, %.10e, %.10e\n", w,AngL(w),Iner(w));
    fflush(fp);
   
    /* Progress bar */ 
    prog = l/(double)maxi;
    printProg(prog);
}
clock_t c1 = clock();
printf("\n------------------------------\n");

/* Bye! */
printf("\n");
printf("Finished calculating %d points in %3.2f sec. \n", maxi, (c1-c0)/1.e6);
printf("\nBYE!\n");
printf("===========================================\n");

return 0;
};

////////////////////////////////////////////////////////////////////////////////
