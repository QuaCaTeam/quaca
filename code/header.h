/* Include the needed libraries */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <complex.h>
#include <string.h>
/* Defining fixed parameters */
#define PI 3.14159265358979323846  // just pi
#define N 3                        // matrix dimension


double complex Ly[3][3] = {
   {0.   , 0., 1.*_Complex_I} ,   /*  initializers for row indexed by 0 */
   {0.   , 0., 0.  } ,   /*  initializers for row indexed by 1 */
   {-1.*_Complex_I, 0., 0.  }     /*  initializers for row indexed by 2 */
};
double complex Ly2[3][3] = {
   {1.   , 0., 0.  } ,   /*  initializers for row indexed by 0 */
   {0.   , 0., 0.  } ,   /*  initializers for row indexed by 1 */
   {0.   , 0., 1.  }     /*  initializers for row indexed by 2 */
};
double  kcut;
double v, za, QF;
double eps0, a0, wa, hbar;
double wp1;
double g1;

// This function multiplies mat1[][] and mat2[][],
// and stores the result in res[][]
void multiply(double complex mat1[N][N], double complex mat2[N][N], double complex res[][N]){
    int i, j, k;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            res[i][j] = 0.;
            for (k = 0; k < N; k++)
                res[i][j] += mat1[i][k]*mat2[k][j];
        }
    }
}

// Yields the transposed of a 3x3 matrix
void dagger(double complex mat1[N][N], double complex mat2[N][N])
{
int i,j;
for (i = 0; i < N; i++)
{
    for (j = 0; j < N; j++)
    {
      mat2[j][i] = conj(mat1[i][j]);
     }
  }
}

void fancyI(double complex mat[N][N], double complex matI[N][N])
{
  int i,j;
  for (i = 0; i < N; i++)
  {
      for (j = 0; j < N; j++)
      {
        matI[i][j] = -0.5*I*( mat[i][j] - conj(mat[j][i]) );
       }
    }
}


// This function yields the trace of a matrix mat
double complex tr(double complex mat[N][N])
{
    int i;
    double complex res;
    res = 0.;
    for (i = 0; i < N; i++)
    {
    res += mat[i][i];
  }
    return res;
}
//==============================================================================
/* INTEGRATION */
// Finite integration int_a^b dx f(x) with a given relative error
double integ ( double my_f() , double a , double b, double relerr)
{
    gsl_function f;
    gsl_integration_cquad_workspace *ws = NULL;
    double res, abserr;
    size_t neval;

    /* Prepare the function. */
    f.function = my_f;
    f.params = NULL;

    /* Initialize the workspace. */
    if ( ( ws = gsl_integration_cquad_workspace_alloc( 200 ) ) == NULL ) {
        printf( "call to gsl_integration_cquad_workspace_alloc failed.\n" );
        abort();
        }

    /* Call the integrator. */
    if ( gsl_integration_cquad( &f, a , b , 0.e-200 , relerr , ws , &res , &abserr , &neval ) != 0 ) {
        printf( "call to gsl_integration_cquad failed.\n" );
        abort();
        }

    /* Free the workspace. */
    gsl_integration_cquad_workspace_free( ws );

    /* Bye. */
    return res;
}

//==============================================================================
/*  FUNCTIONS  */
// Real part of the reflection coefficient r( w , k )
double rR (double w, double k)
{
  double complex rp;
  double complex eps;
  //Z0 = csqrt(cpow(w,2) -cpow(k,2))/w;
  eps = 1. - wp1*wp1/(w*w+w*g1*_Complex_I);
  //Zp = csqrt(epsw*w -cpow(k,2))/(epsw);

rp = (eps  -  1. ) / ( eps + 1. );
return 1.;//creal(rp);
}

// Imaginary part of the reflection coefficient r( w , k )
double rI (double w, double k)
{
double complex rp;
double complex eps;
//Z0 = csqrt(cpow(w,2) -cpow(k,2))/w;
  eps = 1. - wp1*wp1/(w*w+w*g1*_Complex_I);
//Zp = csqrt(epsw*w -cpow(k,2))/(epsw);

rp = (eps  -  1. ) / ( eps + 1. );
return w*g1/(wp1*wp1);//cimag(rp);
}


// Integral over the Green tensor with several options:
void Gint(double complex Gten[N][N], double w, int RorI, int kx, int theta)
{
double Gphi;
double Gsigma[3];
int sig, pphi;
double wrapphi(double phi)
{
  double pre, lim1, lim2, cosp;
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
    return 2*resk/(2*eps0*pow(2*PI,2));
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
  if (theta == 1) {
    if (cosp < 0) {
      if (-w/(v*cosp)<= kcut/(2*za)) {
        lim2 =-w/(v*cosp);
    //    printf("here\n");
      }
    }
  }
  return pre*integ(wrapk, lim1, lim2, 1e-6);
}
// First, we calculate Gsigma
pphi = 0;
// ... the (1,1) component with cos^2(phi) as prefactor
 sig = 0;
 Gsigma[0] = integ( wrapphi, 0., PI, 1e-4);
// ... the (2,2) component with sin^2(phi) as prefactor
  sig = 1;
 Gsigma[1] = integ( wrapphi, 0., PI, 1e-4);
// and the (2,2) component as sin^2 = 1 - cos^2
 Gsigma[2] = Gsigma[0] + Gsigma[1];

// Second, we calculate Gphi
pphi = 1;
 Gphi = integ( wrapphi, 0., PI, 1e-4);
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
void alpha(double complex alp[N][N], double w)
{
double complex a,b,c,d;
double complex alpinv[3][3];
double complex GI[3][3], GR[3][3];

alpinv[0][0] = pow(wa,2)-pow(w,2);
alpinv[1][1] = pow(wa,2)-pow(w,2);
alpinv[2][2] = pow(wa,2)-pow(w,2);

Gint(GI, w, 0, 0, 0);
Gint(GR, w, 1, 0, 0);

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

double IntQF( double w)
{
double complex GIth[3][3], GIk[3][3], GIkth[3][3];
double complex alp[3][3], alpI[3][3], S[3][3],  temp1[3][3], temp2[3][3];
double complex alpdag[3][3];


/* Creating all needed matrices */
alpha(alp,w);
dagger(alp,alpdag);
fancyI(alp, alpI);
Gint(GIth , w, 0, 0, 1);
Gint(GIk  , w, 0, 1, 0);
Gint(GIkth, w, 0, 1, 1);

/* Building the power spectrum S */
multiply(alp,GIth,temp1);
multiply(temp1,alpdag,S);

/* Calculating the trace and return */
multiply(S,GIk,temp1);
multiply(alpI,GIkth,temp2);
//printf("%.10e\n",w );
//printf("%.10e %.10e %.10e\n",creal(tr(temp1)), creal(tr(temp2)) ,creal(tr(temp1)));
return (2*hbar/PI)*creal( -tr(temp1) + tr(temp2) );
}
