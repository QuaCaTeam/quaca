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
double sta = 1E-6 ;                 // start value of the calculation,
double sto = 1E-2 ;                 // final value of the calculation
double spac = pow(sto/sta,1./maxi); // and the respective spacing
double F0 ;                         // normalization constant for QF
double Fanat;                        // analytical approximation (translation)
double Fanar;                        // analytical approximation (rotation)
FILE *fp;                           // output file
int l;                              // dummy index

/* System parameters (input routine is not implemented yet) */
v    = 1E-4;
za   = 5E-9/(1.9732705e-7);  // height of the dipole in 1/eV
eps0 = 1./(4*PI);             // vacuum permittivity
hbar = 1.;                    // reduced Planck's constant
c    = 1.;                    // speed of light
a0   = 6e-9;                  // static polarizability
wa   = 1.3e0;                 // dipole resonance frequency in eV
einf = 3.7;                   // background permittivity
wp1  = 1.30077;               // plasma frequency in eV
wsp1 = wp1/sqrt(1.+einf);     // plasma frequency in eV
g1   = 0.12;                   // damping of the material in eV
kcut = 50.;                   // Integration cut-off of the k-integration
relerr = 1e-2;                // aimed relative error of the integration
recerr = 1e-2;                // increase of relerr per layer of integration
beta   = 1./((3e0)/1.16e4);       // inverse temperature in eV

/* open the file */
fp = fopen("../output/resultsT.dat", "w");
if (fp == NULL) {
   printf("I couldn't open results.dat for writing.\n");
   exit(0);
}


/* Starting calculations */
for (l=0; l<=maxi; ++l){
   printf("progress %3.2f\n",l*100./maxi );
   clock_t c0 = clock();
   /* Point of evaluation */
   v = sta*pow(spac,l);
   /* Performing calculations */
   /* translational contribution */
   transroll = 0;
   absr = 1E-200;
   QFt = integ(IntQF,0.,0.9*wa,relerr, absr);
   QFt = QFt + integ(IntQF,0.5*wa,1.5*wa,relerr, abs(QFt));
   QFt = QFt + integ(IntQF,1.5*wa,100*wa,relerr, abs(QFt));
   /* rolling contribution */
   transroll = 1;
   QFr = integ(IntQF,0.,wa,relerr, absr);
//   QFr = QFr + integ(IntQF,wa,100*wa,relerr, absr);
      /* Calculate normalization constant */
   F0 = -3*hbar*pow(wsp1,5)*a0/(2*PI*eps0*pow(c,4));
   /* Calculating analytical approximation for small velocities */
   Fanat = -63*hbar*a0*a0*pow(g1/(eps0*wp1*wp1),2)*pow(v,3)/(pow(PI,3)*pow(2*za,10));
   Fanar = -45*Fanat/63.;
   Fanat = Fanat - a0*a0*pow(g1/(eps0*wp1*wp1),2)*6*v/(PI*beta*beta*hbar*pow(2*za,8));
   Fanar = Fanar + 3*a0*a0*pow(g1/(eps0*wp1*wp1),2)*v/(PI*beta*beta*hbar*pow(2*za,8));
   /* Print result to the screen */
   printf("v= %.5e\n", v);
   printf("F= %.5e\n", (QFt+QFr)/F0 );
   /* write to the file */
   fprintf(fp, "%.10e, %.10e, %.10e, %.10e, %.10e\n", v, QFt/F0, QFr/F0,Fanat/F0, Fanar/F0);
   /* buffer data for interative writing process */
   fflush(fp);
   clock_t c1 = clock();
   printf("time in sec: %3.2f\n",(c1 - c0) / 1000000. );
   printf ("%s \n", " ");
 }

return 0;
}
////////////////////////////////////////////////////////////////////////////////
