/*! 
 * \file plate.c
 * \brief Functions for quantum friction of a plate.
 * \author M. O.
 *
 * Functions, for more explicit documentation see header.
 */

void fancyI(double complex mat[Ndim][Ndim], double complex matI[Ndim][Ndim])
{
  int i,j;
  for (i = 0; i < Ndim; i++)
  {
      for (j = 0; j < Ndim; j++)
      {
        matI[i][j] = -0.5*_Complex_I*( mat[i][j] - conj(mat[j][i]) );
       }
    }
}

//==============================================================================
/*  FUNCTIONS  */
// Real part of the  reflection coefficient r( w , k )
double rR (double w, double k)
{
  double complex rp;
  double complex epsw/*, Z0, Zp*/;
//  Z0 = csqrt(cpow(w,2) -cpow(k,2))/w;
  epsw = einf*w - wp1*wp1/(w+g1*_Complex_I);
//  Zp = csqrt(epsw*w -cpow(k,2))/(epsw);

//rp = (Z0  -  Zp ) / ( Z0 + Zp );
rp = (epsw - w)/(epsw + w);
return 1.;//creal(rp);
}

// Imaginary part of the reflection coefficient r( w , k )
double rI (double w, double k)
{
  double complex rp;
  double complex epsw/*, Z0, Zp*/;
//  Z0 = csqrt(cpow(w,2) -cpow(k,2))/w;
  epsw = einf*w - wp1*wp1/(w+g1*_Complex_I);
//  Zp = csqrt(epsw*w -cpow(k,2))/(epsw);

//rp = (Z0  -  Zp ) / ( Z0 + Zp );
rp = (epsw - w)/(epsw + w);
return w*g1/(wp1*wp1);//cimag(rp);
}


// Integral over the Green tensor with several options:
void Gint(double complex Gten[Ndim][Ndim], double w, int RorI, int kx, int theta, int T)
{
double Gphi;
double Gsigma[3];
int sig, pphi;
double wrapphi(double phi)
{
  double pre, lim1, lim2, cosp, resphi;
  int caseT=0;
  // Defining the integrand of the k integration
   double wrapk (double k)
   {
   double wpl =  w+cosp*k*v;
   double resk;

   if (RorI == 0) {
     resk = rR(wpl, k)*k*k*exp(-2*za*k);
   }
   if (RorI == 1) {
     resk = rI(wpl, k)*k*k*exp(-2*za*k);
    }
    if (kx == 1) {
     resk = resk*k*cosp;
    }
    if (T == 1) {
      resk = resk/(1.-exp(-beta*hbar*wpl));
    }

    return resk;
  }
  // Performing the k integration
  cosp = cos(phi);
  if (sig == 0) {
  pre = pow(cosp,2);
  }
  if (sig == 1) {
   pre = pow(sin(phi),2);
  }
  if (pphi == 1) {
   pre = cosp;
  }
  lim1 = 0.;
  lim2 = kcut/(2*za);
          caseT = 0;
  if (theta == 1) {
    if (cosp < 0) {
      if (-w/(v*cosp)<= kcut/(2*za)) {
        lim2 =-w/(v*cosp);
        caseT = 1;
    //    printf("here\n");
      }
    }
  }
    resphi = pre*integ(wrapk, lim1, lim2, relerr*recerr*recerr, absr*recerr*recerr);
  if (T == 1 && caseT == 1){
    resphi = resphi + pre*integ(wrapk, lim2, kcut/(2*za), relerr*recerr*recerr, absr*recerr*recerr);
  }

  return resphi/(eps0*pow(2*PI,2));
}
// First, we calculate Gsigma
pphi = 0;
// ... the (1,1) component with cos^2(phi) as prefactor
 sig = 0;
 Gsigma[0] = integ( wrapphi, 0., PI, relerr*recerr, absr*recerr);
// ... the (2,2) component with sin^2(phi) as prefactor
  sig = 1;
 Gsigma[1] = integ( wrapphi, 0., PI, relerr*recerr, absr*recerr);
// and the (2,2) component as sin^2 = 1 - cos^2
 Gsigma[2] = Gsigma[0] + Gsigma[1];

// Second, we calculate Gphi
pphi = 1;
 Gphi = integ( wrapphi, 0., PI, relerr*recerr, absr*recerr);
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
double complex a,b,c,d;
double complex alpinv[3][3];
double complex GI[3][3], GR[3][3];

alpinv[0][0] = wa*wa-w*w;
alpinv[1][1] = alpinv[0][0];
alpinv[2][2] = alpinv[0][0];

Gint(GI, w, 1, 0, 0, 0);
Gint(GR, w, 0, 0, 0, 0);

alpinv[0][0] += -a0*pow(wa,2)*(GR[0][0]+_Complex_I*GI[0][0]);
alpinv[1][1] += -a0*pow(wa,2)*(GR[1][1]+_Complex_I*GI[1][1]);
alpinv[2][2] += -a0*pow(wa,2)*(GR[2][2]+_Complex_I*GI[2][2]);
alpinv[2][0] =  -a0*pow(wa,2)*(GR[2][0]+_Complex_I*GI[2][0]);
alpinv[0][2] = -alpinv[2][0];

  a = alpinv[0][0];
  b = alpinv[1][1];
  c = alpinv[2][2];
  d = alpinv[2][0];

 alp[0][0] = c*a0*pow(wa,2)/(a*c+cpow(d,2));
 alp[0][2] = d*a0*pow(wa,2)/(a*c+cpow(d,2));
 alp[1][1] = a0*pow(wa,2)/b;
 alp[2][0] = -alp[0][2];
 alp[2][2] = a*a0*pow(wa,2)/(a*c+cpow(d,2));
 alp[1][0] = 0.;
 alp[0][1] = 0.;
 alp[1][2] = 0.;
 alp[2][1] = 0.;
}
double anaAngL(double v){
double complex w0[4], sumpol=0.,prodpol=1.;
double G, D;
int i,j;
D  = a0*wa*wa/(4*PI*eps0*pow(2*za,3)) ;
G  = D*2*eps0*g1/(eps0*wp1*wp1);
w0[0] = -_Complex_I*G/2.+csqrt(wa*wa-D-G*G/4.);
w0[1] = -_Complex_I*G/2.-csqrt(wa*wa-D-G*G/4.);
w0[2] = -_Complex_I*G+csqrt(wa*wa-2*D-G*G);
w0[3] = -_Complex_I*G-csqrt(wa*wa-2*D-G*G);

for (i = 0; i < 4; i++)
{
  prodpol = 1.;
  for (j = 0; j < 4; j++)
  {
    if (i != j) {
          prodpol *= (w0[i]-w0[j]);
    }
   }
  sumpol += w0[i]*clog(-w0[i]*2*za/c)/prodpol;
 }

return (3*a0*v*wa*wa*(g1/(eps0*wp1*wp1))/(PI*PI*pow(2*za,4)))*creal(sumpol);
}


double AngL(double w){
  double complex GIth[3][3];
  double complex alp[3][3], S[3][3], temp1[3][3];
  double complex alpdag[3][3], alpI[3][3];
  double Ang;
  int minw=0;

  if (w<0.) {
    w = -w;
    minw = 1;
  }
  /* Creating all needed matrices */
  alpha(alp,w);
  fancyI(alp,alpI);
  dagger(alp,alpdag);
  Gint(GIth , w, 1, 0, 1, 1);

  /* Building the power spectrum S */
  multiply(alp,GIth,temp1);
  multiply(temp1,alpdag,S);

  /* Building the angular momentum and the moment of inertia */
  multiply(S,Ly,temp1);

  /* Returning the angular momentum and the moment of inertia */
  Ang = creal(tr(temp1));

 /* the alpI contribution */
 if (minw==1) {
   multiply(alpI,Ly,temp1);
   Ang = Ang - creal(tr(temp1));
 }

 Ang = (hbar/PI)*w*Ang/(PI*a0*wa*wa);
 if (minw==1){
   w=-w;
 }
 return Ang;
}

double Iner(double w){
  double complex GIth[3][3];
  double complex alp[3][3], S[3][3], temp1[3][3], temp2[3][3];
  double complex alpdag[3][3], alpI[3][3];
  double Ine;
  int minw=0;

  if (w<0.){
    w=-w;
    minw = 1;
  }
  /* Creating all needed matrices */
  alpha(alp,w);
  fancyI(alp,alpI);
  dagger(alp,alpdag);
  Gint(GIth , w, 1, 0, 1, 1);

  /* Building the power spectrum S */
  multiply(alp,GIth,temp1);
  multiply(temp1,alpdag,S);

  /* Building the angular momentum and the moment of inertia */
  multiply(S,Ly2,temp2);

  /* Returning the angular momentum and the moment of inertia */
  Ine = creal(tr(temp2));

 /* the alpI contribution */
 if (minw==1){
    multiply(alpI,Ly2,temp2);
    Ine = Ine - creal(tr(temp2));
    w=-w;
 }


 Ine = (hbar/PI)*Ine/(PI*a0*wa*wa);
 return Ine;
}


double IntQF(double w)
{
double complex GIth[3][3], GIk[3][3], GIkth[3][3];
double complex alp[3][3], alpI[3][3], S[3][3],  temp1[3][3], temp2[3][3];
double complex alpdag[3][3];


/* Creating all needed matrices */
alpha(alp,w);
dagger(alp,alpdag);
fancyI(alp, alpI);
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
//printf("%.10e\n",w );
//printf("%.10e %.10e %.10e\n",creal(tr(temp1)), creal(tr(temp2)) ,creal(tr(temp1)));
return (2*hbar/PI)*creal( -tr(temp1) + tr(temp2) );
}
