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
            else if (strcmp(tmpstr1,"bet1")==0) {
                 bet1 = atof(tmpstr2);
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
               else if (strcmp(tmpstr1,"vF")==0) {
                    vF = atof(tmpstr2);
            }
               else if (strcmp(tmpstr1,"aF")==0) {
                    aF = atof(tmpstr2);
            }
               else if (strcmp(tmpstr1,"Delt")==0) {
                    Delt = atof(tmpstr2);
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
    printf("\na0  : %.5e",a0);

    printf("\n");
    printf("\n// material parameters\n");
    printf("\neinf : %.5e",einf);
    printf("\nwp1  : %.5e in eV",wp1);
    printf("\ng1   : %.5e in eV",g1);
    printf("\nbet1 : %.5e      ",bet1);
    printf("\nT    : %.5e in K",T);

    printf("\n");
    printf("\n// numerical specifications\n");
    printf("\nkcut   : %.5e",kcut);
    printf("\nrelerr : %.5e",relerr);
    printf("\nrecerr : %.5e",recerr);
    printf("\nabserr : %.5e",abserr);

    printf("\n");
    printf("\n// graphene specifications\n");
    printf("\nvF   : %.5e",vF);
    printf("\naF   : %.5e",aF);
    printf("\nDelt : %.5e",Delt);

    printf("\n");
    printf("\n===========================================\n");
// transforming the remaining SI units to natural units
// or other convenient forms
za   = za/(1.9732705e-7);
beta = 1./(T/1.16e4);
wsp1 = wp1/sqrt(1.+einf);
OmD  = 2*Delt;
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
 if (w<0) {
   r[0] = conj(r[0]);
   r[1] = conj(r[1]);
 }
}

void reflhydro(double complex r[2] , double w, double complex kap)
{
 double complex epsw, kapD, kapL, k2, kapDLw;
 double wabs;
 // Calculating the dielectric function (times w due to numerical convenience)
 wabs = fabs(w);
 k2 = kap*kap+w*w;
 epsw = einf*wabs - wp1*wp1/(wabs+g1*_Complex_I);

 // Calculating the wavevector in z direction...
 kapD   = csqrt(kap*kap-(epsw-wabs)*wabs);
 // ... and fixing the signs
 kapD   = fabs(creal(kapD)) - fabs(cimag(kapD))*_Complex_I;

 kapL   = csqrt(kap*kap+wabs*wabs+(wp1*wp1-wabs*(wabs+_Complex_I*g1))/(bet1*bet1*vF*vF) );
 kapL   = fabs(creal(kapL)) - fabs(cimag(kapL))*_Complex_I;
 kapDLw = kapD*wabs + (epsw-wabs)*k2/kapL;
// Calculating r^s
 r[1]   =  -(kapD  -  kap ) / ( kapD + kap );
// Calculating r^p
 r[0]   = (epsw*kap  -  kapDLw ) / ( epsw*kap + kapDLw );
 if (w<0) {
   r[0] = conj(r[0]);
   r[1] = conj(r[1]);
 }
}

void reflGraph(double complex r[2] , double w, double complex kap)
{
 double complex kapF, wabs;
 double complex kapFOm, PP00, PPtr, Phi;
 double k2, kapF2, kap2;
 // Calculating the dielectric function (times w due to numerical convenience)
 wabs = fabs(w);
 kap2 = creal(kap*kap);
 k2 = creal(kap*kap) + wabs*wabs ;
 // Calculating the wavevector in z direction...
 kapF   = csqrt(vF*vF*k2-wabs*wabs);
 // ... and fixing the signs
 kapF  = fabs(creal(kapF)) - fabs(cimag(kapF))*_Complex_I;
 kapF2 = creal(kapF*kapF);

 kapFOm = kapF/OmD;

 Phi = 2*OmD*(1.+( kapFOm - 1./kapFOm )*catan(kapFOm)  );
 PP00 = Phi*aF*k2/kapF2;
 PPtr = Phi*aF*(kap2+kapF2)/kapF2;

// Calculating r^s
 r[1]   =  - (k2*PPtr-kap2*PP00)/(k2*(PPtr+2*kap)-kap2*PP00);
// Calculating r^p
 r[0]   = kap*PP00/(kap*PP00+2*k2);
 if (w<0) {
   r[0] = conj(r[0]);
   r[1] = conj(r[1]);
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
  double limI1, limgraph;
  int caseT=0;
  // Defining the integrand of the k integration
   double wrapkap (double kap)
   {
   double wpl;
   double cosp2 = cosp*cosp;
   double sinp2 = sinp*sinp;
   double resk=0., k;
   double complex kapc, r[2], resc, rppre, rspre;
   double v2 = v*v ;
   double kapc2;
   double complex prefac;
   if (kap < 0) {
     kapc = _Complex_I*kap + 0.;
   }
   else if (kap >= 0) {
     kapc = kap + _Complex_I*0.;
   }
   kapc2 = creal(kapc*kapc);
   k   = ((sqrt(kapc2*(1.-v2*cosp2) +w*w) +v*w*cosp)/(1.-v2*cosp2) );
   wpl =  w+cosp*k*v;
   refl(r,wpl,kapc);

   prefac = cexp(-2*za*kapc)/(1.-cosp*v*wpl/k);

   rppre = r[0]*prefac;
   rspre = r[1]*prefac;

   // Here we are calculating the phi(k,wpl)
   if (pphi == 1) {
    resc = rppre*k*kapc*cosp;
   }
   // Here we are calculating the different sigma(k,wpl) components
   if (pphi == 0) {
    rspre = rspre*wpl*wpl;
    if (sig == 0) {
      resc = (rppre*cosp2*kapc2 + rspre*sinp2);
    }
    if (sig == 1) {
      resc = (rppre*sinp2*kapc2 + rspre*cosp2);
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
//    if (T == 1) {
//      resk = resk/(1E0-exp(-beta*wpl));
//    }
    return resk;
  }
  // Performing the k integration
  cosp = cos(phi);
  sinp = sin(phi);
//  limgraph = sqrt(OmD*OmD +(1.-vF*vF)*k*k);
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



//  if (T == 1 && caseT == 1){
//    resphi += integ(wrapkap, lim2,kcut/(2*za), relerr*recerr*recerr, abserr*recerr*recerr);
//  }

  return resphi/PI;
}
// First, we calculate Gsigma
pphi = 0;
// ... the (1,1) component with cos^2(phi) as prefactor
  sig = 0;
  Gsigma[0] = integ( wrapphi, 0., PI, relerr*recerr, abserr*recerr);
// ... the (2,2) component with sin^2(phi) as prefactor
  sig = 1;
  Gsigma[1] = integ( wrapphi, 0., PI, relerr*recerr, abserr*recerr);
// and the (2,2) component as sin^2 = 1 - cos^2
  sig = 2;
  Gsigma[2] = integ( wrapphi, 0., PI, relerr*recerr, abserr*recerr);

// Second, we calculate Gphi
pphi = 1;
 Gphi = integ( wrapphi, 0., PI, relerr*recerr, abserr*recerr);
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

alpinv00 = (wa*wa-w*w)/(a0*wa*wa);
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


double IntQF( double w)
{
double complex GIth[3][3], GIk[3][3], GIkth[3][3];
double complex alp[3][3], alpI[3][3], S[3][3],  temp1[3][3], temp2[3][3];
double complex alpdag[3][3];


/* Creating all needed matrices */
alpha(alp,w);
dagger(alp,alpdag);
fancy(alp, alpI, -1);
Gint(GIth , w, 1, 0, 1, 1);
Gint(GIk  , w, 1, 1, 0, 0);
Gint(GIkth, w, 1, 1, 1, 1);

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

return (2./PI)*creal( -tr(temp1) + tr(temp2) );
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
