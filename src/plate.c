/*!
 * \file plate.c
 * \brief Functions for quantum friction of a plate.
 * \author M. O.
 *
 * Functions, for more explicit documentation see header.
 */

/* --------- */
/* Libraries */
/* --------- */

#include "h/qfhelp.h"
#include "h/plate.h"

/* --------- */
/* Functions */
/* --------- */

void input()
{
FILE * fr = fopen("../data/input/paratest.dat", "rt");


if(fr == NULL){printf("file %s not found","../data/input/paratest.dat");}

char tmpstr1[16];
char tmpstr2[16];

double T;


    char tempbuff[100];

    while(!feof(fr))
    {
         if (fgets(tempbuff,100,fr)) {
            sscanf(tempbuff, "%15s : %15s", tmpstr1, tmpstr2);

            if (strcmp(tmpstr1,"wp1")==0) {
                 wp1 = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"g1")==0) {
                 g1 = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"einf")==0) {
                 einf = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"v")==0) {
                 v = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"za")==0) {
                 za = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"T")==0) {
                 T = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"a0")==0) {
                 a0 = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"kcut")==0) {
                 kcut = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"relerr")==0) {
                 relerr = atof(tmpstr2);
            }
               else if (strcmp(tmpstr1,"abserr")==0) {
                 abserr = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"recerr")==0) {
                 recerr = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"wa")==0) {
                 wa = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"gamMu")==0) {
                 gamMu = atof(tmpstr2);
            }
	    else if (strcmp(tmpstr1,"muquest")==0) {
                 muquest = atoi(tmpstr2);
            }
               else if (strcmp(tmpstr1,"vF")==0) {
                    vF = atof(tmpstr2);
            }
               else if (strcmp(tmpstr1,"aF")==0) {
                    aF = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1,"//")==0) {
                 /* skip */
            }
            else{
                printf("Unrecongized parameter : \"%s\"\n", tmpstr1);
            }
         }
    }
    fclose(fr);

// Recall the read parameters and print to screen
    printf("\n===========================================\n");
    printf("\nINPUT PARAMETERS:");
    printf("\n`````````````````\n");
    printf("\n");
    printf("// atomic constants\n");
    printf("\nv   : %.5e in c",v);
    printf("\nza  : %.5e in nm",za);
    printf("\nwa  : %.5e in eV",wa);
    printf("\nga  : %.5e in eV",gamMu);
    printf("\nmu? : %i",muquest);
    printf("\na0  : %.5e",a0);

    printf("\n");
    printf("\n// material parameters\n");
    printf("\neinf : %.5e",einf);
    printf("\nwp1  : %.5e in eV",wp1);
    printf("\ng1   : %.5e in eV",g1);
    printf("\nvF   : %.5e",vF);
    printf("\naF   : %.5e",aF);
    printf("\nT    : %.5e in K",T);

    printf("\n");
    printf("\n// numerical specifications\n");
    printf("\nkcut   : %.5e",kcut);
    printf("\nrelerr : %.5e",relerr);
    printf("\nrecerr : %.5e",recerr);
    printf("\nabserr : %.5e",abserr);


    printf("\n");
    printf("\n===========================================\n");
// transforming the remaining SI units to natural units
// or other convenient forms
za   = za/(1.9732705e-7);
beta = 1./(T/1.16e4);
wsp1 = wp1/sqrt(1.+einf);
//

}


void refl(double complex r[2] , double w, double complex kap)
{
 double complex epsw, kaplc, wabs;
 // Calculating the dielectric function (times w due to numerical convenience)
 wabs = fabs(w);
 epsw = einf*wabs - wp1*wp1/(wabs+g1*_Complex_I);

 // Calculating the wavevector in z direction...
 kaplc   = csqrt(kap*kap-(epsw-wabs)*wabs);
 // ... and fixing the signs
 kaplc   = fabs(creal(kaplc)) - fabs(cimag(kaplc))*_Complex_I;

// Calculating r^s
 r[1]   =  -(kaplc  -  kap ) / ( kaplc + kap );
// Calculating r^p
 r[0]   = (epsw*kap  -  wabs*kaplc ) / ( epsw*kap + wabs*kaplc );
//  r[0]   = (epsw  -  wabs ) / ( epsw + wabs );
 if (w<0) {
   r[0] = conj(r[0]);
   r[1] = conj(r[1]);
 }
}

void reflhydro(double complex r[2] , double w, double complex kap)
{
 double complex epsw, kapD, kapL, k2, kapDLw, bet2, kap2;
 double wabs;
 // Calculating the dielectric function (times w due to numerical convenience)
 wabs = fabs(w);
 bet2 =( (3./5.)*wabs+(1./3.)*g1*_Complex_I )/(wabs+g1*_Complex_I);
 kap2 = kap*kap;
 k2 = kap2+w*w;

 epsw = einf*wabs - wp1*wp1/(wabs+g1*_Complex_I);

 // Calculating the wavevector in z direction...
 kapD   = csqrt(kap2-(epsw-wabs)*wabs);
 // ... and fixing the signs
 kapD   = fabs(creal(kapD)) - fabs(cimag(kapD))*_Complex_I;

 kapL   = csqrt(kap2+wabs*wabs+(wp1*wp1-wabs*(wabs+_Complex_I*g1))/bet2 );
 kapL   = fabs(creal(kapL)) - fabs(cimag(kapL))*_Complex_I;
 kapDLw = kapD*wabs + (epsw-wabs)*k2/kapL;
// Calculating r^s
 r[1]   =  -(kapD  -  kap ) / ( kapD + kap );
// Calculating r^p
 r[0]   = (epsw*kap  -  kapDLw ) / ( epsw*kap + kapDLw );
 if ( w < 0 && omsgn == 0) {
   // For the Omega calcualtion, we need a integration of the Green tensor with
   // an additional sgn( kx*v + w ). This is equivalent with calculating
   // r_I( abs(w) ).
   r[0] = conj(r[0]);
 }
}


// Integral over the Green tensor with several options:
void Gint(double complex Gten[Ndim][Ndim], double w, int RorI, int kx, int theta, int T)
{
double Gphi;
double Gsigma[3];
int sig, pphi;
double wrapphi(double phi)
{
  double lim2, cosp, sinp, resphi;
  double limI1;
  int caseT=0;

  // Defining the integrand of the k integration
   double wrapkap (double kap)
   {
   double wpl;
   double cosp2 = cosp*cosp;
   double sinp2 = sinp*sinp;
   double resk=0.;
   double complex kapc=0., re[2], resc=0., rppre, rspre;
   double v2 = v*v, k ;
   double kap2 = kap*fabs(kap);
   double complex prefac;
   if (kap >= 0E0) {
     kapc = kap + I*0E0;
   }else {
     kapc = 0E0 + I*kap;
   }
   k = (sqrt(kap2*(1E0-v2*cosp2) + w*w) + v*w*cosp)/(1E0-v2*cosp2);
   wpl =  w+cosp*k*v;
   re[0] = 0E0;
   re[1] = 0E0;
   refl(re,wpl,kapc);
   prefac = cexp(-2*za*kapc)/(1E0-cosp*v*wpl/k);

   rppre = re[0]*prefac;
   rspre = re[1]*prefac;

   // Here we are calculating the phi(k,wpl)
   if (pphi == 1) {
    resc = rppre*cosp*k*kapc;
   }
   // Here we are calculating the different sigma(k,wpl) components
   if (pphi == 0) {
      rspre = rspre*wpl*wpl;
    if (sig == 0) {
      resc = (rppre*cosp2*kap2 + rspre*sinp2);
    }
    if (sig == 1) {
      resc = (rppre*sinp2*kap2 + rspre*cosp2);
    }
    if (sig == 2) {
      resc = rppre*k*k;
    }
   }

   if (RorI == 0) {
     resk = creal(resc);
   }
   if (RorI == 1) {
     resk = cimag(resc);
    }
    if (kx == 1) {
     resk = resk*k*cosp;
    }
    if (T == 1) {
      resk = resk*(1E0/(1E0-exp(-beta*wpl)) - 1E0/(1E0-exp(-beta*w)) );
    }
    return resk;
  }
  // Performing the kappa integration
  cosp = cos(phi);
  sinp = sin(phi);

  limI1 = -w;
  lim2 = kcut/(2*za);
          caseT = 0;
  if (theta == 1) {
    if (cosp < 0) {
      if (-w/(v*cosp)< kcut/(2*za)) {
        lim2 =-w/(v*cosp);
        caseT = 1;
      }
    }
  }

    resphi = integ(wrapkap,limI1, lim2, relerr*recerr*recerr, abserr*recerr*recerr);

  if (T == 1 && caseT == 1){
    resphi += integ(wrapkap, lim2,kcut/(2*za), relerr*recerr*recerr, abserr*recerr*recerr);
  }
  return resphi/PI;
}
// First, we calculate Gsigma
pphi = 0;
// ... the (1,1) component with cos^2(phi) as prefactor
  sig = 0;
  Gsigma[0] = integ( wrapphi, 0., PI*0.5, relerr*recerr, abserr*recerr);
  Gsigma[0] += integ( wrapphi, PI*0.5, PI, relerr*recerr, abserr*recerr);
// ... the (2,2) component with sin^2(phi) as prefactor
  sig = 1;
  Gsigma[1] = integ( wrapphi, 0., PI*0.5, relerr*recerr, abserr*recerr);
  Gsigma[1] += integ( wrapphi, PI*0.5, PI, relerr*recerr, abserr*recerr);
// and the (2,2) component as sin^2 = 1 - cos^2
  sig = 2;
  Gsigma[2] = integ( wrapphi, 0., PI*0.5, relerr*recerr, abserr*recerr);
  Gsigma[2] += integ( wrapphi, PI*0.5, PI, relerr*recerr, abserr*recerr);

// Second, we calculate Gphi
pphi = 1;
 Gphi = integ( wrapphi, 0., PI*0.5, relerr*recerr, abserr*recerr);
 Gphi += integ( wrapphi, PI*0.5, PI, relerr*recerr, abserr*recerr);
// Now we can assemble the full Green tensor
Gten[0][0] = Gsigma[0];
Gten[1][1] = Gsigma[1];
Gten[2][2] = Gsigma[2];
Gten[0][2] = -_Complex_I*Gphi;
Gten[2][0] = -Gten[0][2];
Gten[0][1] = 0.;
Gten[1][0] = 0.;
Gten[2][1] = 0.;
Gten[1][2] = 0.;
}


// The polarizability
void alpha(double complex alp[Ndim][Ndim], double w)
{
double complex alpinv00,alpinv11,alpinv22,alpinv20, det;
double complex GI[3][3], GR[3][3];

alpinv00 =  wa*wa - w*w ;
 if (muquest == 1) {
	alpinv00 = alpinv00 - w*I*mu(w);
 }
alpinv00 = alpinv00/(a0*wa*wa);

alpinv11 = alpinv00;
alpinv22 = alpinv00;

Gint(GI, w, 1, 0, 0, 0);
Gint(GR, w, 0, 0, 0, 0);

alpinv00 = alpinv00 - ( GR[0][0] + I * GI[0][0] );
alpinv11 = alpinv11 - ( GR[1][1] + I * GI[1][1] );
alpinv22 = alpinv22 - ( GR[2][2] + I * GI[2][2] );
alpinv20 = -( GR[2][0] + I * GI[2][0] );

det = alpinv20*alpinv20 + alpinv00*alpinv22;
alp[0][0] = alpinv22/det;
alp[1][1] = (double complex)1./alpinv11;
alp[2][2] = alpinv00/det;
alp[0][2] = alpinv20/det;
alp[2][0] = -alp[0][2];
alp[1][0] =(double complex) 0.;
alp[0][1] =(double complex) 0.;
alp[2][1] =(double complex) 0.;
alp[1][2] =(double complex) 0.;
}

double complex mu( double w)
{
double complex muresC;
muresC = gamMu +I*0E0;

return muresC;
}	

double IntQF( double w)
{
double complex GIth[3][3], GIk[3][3], GIkth[3][3];
double complex alp[3][3], alpI[3][3], S[3][3],  temp1[3][3], temp2[3][3];
double complex alpdag[3][3];
double mur = w*creal(mu(w))/( a0*wa*wa*(1E0-exp(-beta*w)) );

/* Creating all needed matrices */
alpha(alp,w);
dagger(alp,alpdag);
fancy(alp, alpI, -1);
Gint(GIth , w, 1, 0, 1, 1);
Gint(GIk  , w, 1, 1, 0, 0);
Gint(GIkth, w, 1, 1, 1, 1);
if ( muquest == 1 ){
	GIth[0][0] += mur;
	GIth[1][1] += mur;
	GIth[2][2] += mur;
}

/* Building the power spectrum S */
multiply(alp,GIth,temp1);
multiply(temp1,alpdag,S);

/* Calculating the trace and return */
/* We split the calcualtion in translational and rolling friction */

if (transroll == 0){
GIk[2][0] = 0.;
GIk[0][2] = 0.;
GIkth[2][0] = 0.;
GIkth[0][2] = 0.;
}
if (transroll == 1){
GIk[0][0] = 0.;
GIk[1][1] = 0.;
GIk[2][2] = 0.;
GIkth[0][0] = 0.;
GIkth[1][1] = 0.;
GIkth[2][2] = 0.;
}

multiply(S,GIk,temp1);
multiply(alpI,GIkth,temp2);
printf("\n tr(S)=%.5e",creal(tr(temp1)));
return (2E0/PI)*creal( -tr(temp1)  + tr(temp2) );
}


double Omega( double w)
{

}



double IntQFfree( double w)
{
double complex GIkth[3][3], GIkthT[3][3];
double complex alp[3][3], alpI[3][3],  temp1[3][3], temp2[3][3];
/* Creating all needed matrices */
alpha(alp,w);
fancy(alp, alpI,-1);
Gint(GIkth , w, 1, 1, 1, 0);
Gint(GIkthT, w, 1, 1, 1, 1);

/* Calculating the trace and return */
/* We split the calcualtion in translational and rolling friction */
if (transroll == 0){
GIkth[2][0] = 0.;
GIkth[0][2] = 0.;
GIkthT[2][0] = 0.;
GIkthT[0][2] = 0.;
}
if (transroll == 1){
GIkth[0][0] = 0.;
GIkth[1][1] = 0.;
GIkth[2][2] = 0.;
GIkthT[0][0] = 0.;
GIkthT[1][1] = 0.;
GIkthT[2][2] = 0.;
}

multiply(alpI,GIkthT,temp1);
multiply(alpI,GIkth,temp2);
return (2./PI)*creal( tr(temp1) - tr(temp2) );
}
