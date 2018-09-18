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
int maxi=100;                        // plot points
double sta = 1E-9/(1.9732705e-7) ;                 // start value of the calculation,
double sto = 100E-9/(1.9732705e-7) ;                 // final value of the calculation
double spac = pow(sto/sta,1./maxi);      // and the respective spacing
FILE *fp;                           // output file
int l;                              // dummy index
double resAngL, resIner;
/* System parameters (input routine is not implemented yet) */
v    = 1E-4;                   // velocity in c
za   = 20E-9/(1.9732705e-7);   // height of the dipole in 1/eV
eps0 = 1./(4*PI);             // vacuum permittivity
hbar = 1.;                    // reduced Planck's constant
c    = 1.;                    // speed of light
a0   = 6e-9;                  // static polarizability
wa   = 1.3e0;                 // dipole resonance frequency in eV
einf = 1.;                   // background permittivity
wp1  = 9.;               // plasma frequency in eV
wsp1 = wp1/sqrt(1.+einf);     // plasma frequency in eV
g1   = 0.1;                   // damping of the material in eV
kcut = 70.;                   // Integration cut-off of the k-integration
relerr = 1e-4;                // aimed relative error of the integration
recerr = 1e-2;                // increase of relerr per layer of integration
beta   = 1./((1e-9)/1.16e4);       // temperature in eV
delta  = a0*wa*wa/(4*PI*eps0*pow(2*za,3));

/* open the file */
fp = fopen("../output/resultsOm.dat", "w");
if (fp == NULL) {
   printf("I couldn't open results.dat for writing.\n");
   exit(0);
}


/* Starting calculations */
for (l=0; l<=maxi; ++l){
      clock_t c0 = clock();
   printf("progress %3.2f\n",l*100./maxi );
   /* Point of evaluation */
   za = sta*pow(spac,l);
   /* Print result to the screen */
   printf("za= %.10e, %.10e, %.10e\n", za, sqrt(wa*wa-2*delta), sqrt(wa*wa-delta));
   /* write to the file */
   resAngL = integ(AngL,0,sqrt(wa*wa-2*delta),relerr,absr);
    printf("L1= %.10e, %.10e\n",resAngL,anaAngL(v));
   resAngL = resAngL+integ(AngL,sqrt(wa*wa-2*delta),sqrt(wa*wa-delta),relerr,absr);
    printf("L2= %.10e, %.10e\n",resAngL,anaAngL(v));
   resAngL = resAngL+ integ(AngL,sqrt(wa*wa-delta),wa,relerr,absr);
    printf("L3= %.10e, %.10e\n",resAngL,anaAngL(v));
   resAngL = resAngL+ integ(AngL,wa,wsp1,relerr,absr);
    printf("L4= %.10e, %.10e\n",resAngL,anaAngL(v));
   resAngL = resAngL+ integ(AngL,wsp1,wsp1*100,relerr,absr);
    printf("Lf= %.10e, %.10e\n",resAngL,anaAngL(v));
//   resIner = integ(Iner,0.,wa-2*delta,relerr,absr);
//   resIner += integ(Iner,wa-2*delta,wa-delta,relerr,absr);
//   resIner += integ(Iner,wa-delta,wsp1,relerr,absr);
//   resIner += integ(Iner,wsp1,wsp1*100,relerr,absr);
   fprintf(fp, "%.10e, %.10e, %.10e, %.10e\n", za,resAngL,anaAngL(v),resIner);

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
