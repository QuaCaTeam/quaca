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
#include "header.h"
/******************************************************************************/


/******************************************************************************/
// Main program
int main () {

/* Plot parameters */
int maxi=20;                        // plot points
double sta = 1E-5 ;                 // start value of the calculation,
double sto = 1E-3 ;                 // final value of the calculation
double spac = pow(sto/sta,1./maxi); // and the respective spacing
FILE *fp;                           // output file
int l;                              // dummy index

/* System parameters (input routine is not implemented yet) */
za   = 10E-9/(1.9732705e-7);  // height of the dipole in 1/eV
eps0 = 1./(4*PI);             // vacuum permeabillity
hbar = 1.;                    // reduced Planck's constant
a0   = 1e-9;                  // static polarizability
wa   = 1.3e0;                 // dipole resonance frequency in eV
wp1  = 9.;                    // plasma frequency in eV
g1   = 0.1;                   // damping of the material in eV
kcut = 20.;                   // Integration cut-off of the k-integration

/* open the file */
fp = fopen("../output/results.dat", "w");
if (fp == NULL) {
   printf("I couldn't open results.dat for writing.\n");
   exit(0);
}


/* Starting calculations */
for (l=0; l<=maxi; ++l){

   /* Point of evaluation */
   v = sta*pow(spac,l);
   /* Performing calculations */
   QF = integ(IntQF,0,wp1,1e-4);
   /* Print result to the screen */
   printf("%2.20f, %2.20f\n", v, QF);
   /* write to the file */
   fprintf(fp, "%2.20f, %2.20f\n", v, QF);
   /* buffer data for interative writing process */
   fflush(fp);
 }

return 0;
}
////////////////////////////////////////////////////////////////////////////////
