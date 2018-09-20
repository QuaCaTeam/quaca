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


/******************************************************************************/
// Main program
int main () {

/* Plot parameters */
int maxi=1000;                        // plot points
double sta = 1.29999953 ;                 // start value of the calculation,
double sto = 1.299999533;                 // final value of the calculation
double spac =(sto-sta)/maxi;      // and the respective spacing
FILE *fp;                           // output file
int l;                              // dummy index
double w;
/* System parameters (input routine is not implemented yet) */
v    = 1E-3;                   // velocity in c
za   = 20E-9/(1.9732705e-7);   // height of the dipole in 1/eV
eps0 = 1./(4*PI);             // vacuum permittivity
hbar = 1.;                    // reduced Planck's constant
c    = 1.;                    // speed of light
a0   = 6e-9;                  // static polarizability
wa   = 1.3e0;                 // dipole resonance frequency in eV
einf = 1.;                   // background permittivity
wp1  = 12.;               // plasma frequency in eV
wsp1 = wp1/sqrt(1.+einf);     // plasma frequency in eV
g1   = 0.1;                   // damping of the material in eV
kcut = 100.;                   // Integration cut-off of the k-integration
relerr = 1e-8;                // aimed relative error of the integration
recerr = 1e-3;                // increase of relerr per layer of integration
beta   = 1./((1e-3)/1.16e4);       // temperature in eV
delta  = a0*wa*wa/(4*PI*eps0*pow(2*za,3));

/* open the file */
fp = fopen("../output/resultsOmInt.dat", "w");
if (fp == NULL) {
   printf("I couldn't open results.dat for writing.\n");
   exit(0);
}


/* Starting calculations */
for (l=0; l<=maxi; ++l){
      clock_t c0 = clock();
   printf("progress %3.2f\n",l*100./maxi );
   /* Point of evaluation */
   w = sta+spac*l;
   /* Print result to the screen */
   printf("w= %.10e\n", w);
   /* write to the file */
   fprintf(fp, "%.10e, %.10e, %.10e\n", w,AngL(w),Iner(w));

  // fprintf(fp, "%.10e, %.10e, %.10e\n", w,AngL(w),Iner(w));
   /* buffer data for interative writing process */
   fflush(fp);
   clock_t c1 = clock();
   printf("time in sec: %3.2f\n",(c1 - c0) / 1000000. );
   printf ("%s \n", " ");
 }

return 0;
}
////////////////////////////////////////////////////////////////////////////////
