/* --------- */
/* LIBRARIES */
/* --------- */

#include "qfhelp.h"
#include "plate.h"

/* --------- */
/* FUNCTIONS */
/* --------- */

// input routine
void inputPlate(char file[], unsigned int verbose) {
    FILE * fr = fopen(file, "rt");

    if (fr == NULL) {
        printf("file %s not found\n", file);
        exit(0);
    }


    /* dummies */
    char s[100];
    char line[200];
    double i;

    /* read lines till end of stream */
    while(!feof(fr)) {
        /* if line is read */ 
        if (fgets(line, 200, fr)) {

            /* if string and number immediately after are read */    
            if ( sscanf(line, "%s %lf", s, &i) == 2) {
                if (strcmp(s,"wp1")==0) {
                    inputparams.wp1 = i;
                } else if (strcmp(s,"g1")==0) {
                    inputparams.g1 = i;
                } else if (strcmp(s,"einf")==0) {
                    inputparams.einf = i;
                } else if (strcmp(s,"v")==0) {
                    inputparams.v = i;
                } else if (strcmp(s,"za")==0) {
                    inputparams.za = i;
                } else if (strcmp(s,"T")==0) {
                    inputparams.beta = 1E0/(i*kB);
                } else if (strcmp(s,"a0")==0) {
                    inputparams.a0 = i;
                } else if (strcmp(s,"kcut")==0) {
                    inputparams.kcut = i;
                } else if (strcmp(s,"relerr")==0) {
                    inputparams.relerr = i;
                } else if (strcmp(s,"abserr")==0) {
                    inputparams.abserr = i;
                } else if (strcmp(s,"recerr")==0) {
                    inputparams.recerr = i;
                } else if (strcmp(s,"wa")==0) {
                    inputparams.wa = i;
                } else if (strcmp(s,"vF")==0) {
                    inputparams.vF = i;
                } else if (strcmp(s,"aF")==0) {
                    inputparams.aF = i;
                } else if (strcmp(s,"gamMu")==0) {
                    inputparams.gamMu = i;
                } else if (strcmp(s,"muquest")==0) {
                    inputparams.muquest = (int) i;
                } else if (strcmp(s,"start")==0) {
                    inputparams.start = i;
                } else if (strcmp(s,"stop")==0) {
                    inputparams.stop = i;
                } else if (strcmp(s,"steps")==0) {
                    inputparams.steps = (int) i;
                } else {
                    printf("Unrecongized parameter : \"%s\"\n", s);
                }
            }

            if ( sscanf(line, "runvar %s", s) == 1 ) {
                strcpy(inputparams.runvar, s); 
            } 

            if ( sscanf(line, "scale %s", s) == 1 ) {
                strcpy(inputparams.scale, s); 
            } 
        }
    }


    // recall the read parameters
    if (verbose == 1) {
        printf("INPUT PARAMETERS:\n");
        printf("// atomic constants\n");
        printf("v   = %.5e\n",inputparams.v);
        printf("za  = %.5e\n",inputparams.za);
        printf("wa  = %.5e\n",inputparams.wa);
        printf("a0  = %.5e\n",inputparams.a0);
        printf("\n// material parameters\n");
        printf("einf = %.5e\n",inputparams.einf);
        printf("wp1  = %.5e\n",inputparams.wp1);
        printf("g1   = %.5e\n",inputparams.g1);
        printf("vF   = %.5e\n",inputparams.vF);
        printf("aF   = %.5e\n",inputparams.aF);
        printf("T    = %.5e\n",1.16e4/inputparams.beta);
        printf("\n// numerical specifications\n");
        printf("kcut   = %.5e\n",inputparams.kcut);
        printf("relerr = %.5e\n",inputparams.relerr);
        printf("recerr = %.5e\n",inputparams.recerr);
        printf("abserr = %.5e\n",inputparams.abserr);
        printf("\n\nCalculating %s from %.5e to %.5e with %d points on %s scale.\n", inputparams.runvar, inputparams.start, inputparams.stop, inputparams.steps, inputparams.scale);
    }

    // transforming the remaining SI units to natural units
    // or other convenient forms
    inputparams.za   = inputparams.za/(hbar*cvac);
    inputparams.wsp1 = inputparams.wp1/sqrt(1.+inputparams.einf);
    inputparams.eps0= 1./(4*PI);
}

// reflection coefficient
void refl(double complex r[2], double w, double complex kap, void * p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double einf = (params->einf);
    const double wp1 = (params->wp1);
    const double g1 = (params->g1);

    double complex epsw, kaplc;
    double wabs;
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

// Integral over the Green tensor with several options:
void Gint(double complex Gten[Ndim][Ndim], double w, void * p, unsigned int RorI, unsigned int kx, unsigned int theta, unsigned int T) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double v = (params->v);
    const double za = (params->za);
    const double beta = (params->beta);
    const double relerr = (params->relerr);
    const double recerr = (params->recerr);
    const double abserr = (params->abserr);
    const double kcut = (params->kcut);

    double Gphi;
    double Gsigma[3];
    double complex re[2];
    register unsigned int sig, pphi;

    double lim2, limI1;
    limI1 = -w;
    lim2 = kcut/(2.*za);

    double v2 = v*v;

    double wrapphi(double phi) {
        double cosp, sinp, resphi;
        register unsigned int caseT=0;
        
        // Performing the kappa integration
        cosp = cos(phi);
        sinp = sin(phi);

        double cosp2 = cosp*cosp;
        double sinp2 = sinp*sinp;

        // Defining the integrand of the k integration
        double wrapkap (double kap) {
            double wpl;
            double resk=0.;
            double complex kapc=0., resc=0., rppre, rspre;
            double k;
            double kap2 = kap*fabs(kap);
            double complex prefac;
            if (kap >= 0E0) {
                kapc = kap + I*0E0;
            } else {
                kapc = 0E0 + I*kap;
            }
            k = (sqrt(kap2*(1E0-v2*cosp2) + w*w) + v*w*cosp)/(1E0-v2*cosp2);
            wpl =  w+cosp*k*v;
            refl(re,wpl,kapc,params);
            prefac = cexp(-2*za*kapc)/(1E0-cosp*v*wpl/k);

            rppre = re[0]*prefac;
            rspre = re[1]*prefac;

            // Here we are calculating the phi(k,wpl)
            switch(pphi){
                case 0:
                    rspre = rspre*wpl*wpl;

                    switch (sig) {
                        case 0:
                            resc = (rppre*cosp2*kap2 + rspre*sinp2);
                            break;
                        case 1:
                            resc = (rppre*sinp2*kap2 + rspre*cosp2);
                            break;
                        case 2:
                            resc = rppre*k*k;
                            break;
                    }
                    break;
                case 1:
                    resc = rppre*cosp*k*kapc;
                    break;
            }

            switch(RorI) {
                case 0:
                    resk = creal(resc);
                    break;
                case 1:
                    resk = cimag(resc);
                    break;
            }

            if (kx == 1) {
                resk = resk*k*cosp;
            }

            if (T == 1) {

		    resk = resk*(
			1E0/(1E0-exp(-beta*wpl)) - 1E0/(1E0-exp(-beta*w)) 
			     	);
	    }     
            return resk;
        }


        caseT = 0;
        if (theta == 1 && cosp < 0 && (-w/(v*cosp)< kcut/(2*za)) ) {
            lim2 =-w/(v*cosp);
            caseT = 1;
        }

        resphi = cquad(wrapkap, params, limI1, lim2, relerr*recerr*recerr, abserr*recerr*recerr);
	
	if (T == 1 && caseT == 1){
        
	    	resphi += cquad(wrapkap, params, lim2, kcut/(2.*za), relerr*recerr*recerr, abserr*recerr*recerr);
       
       	}
        return resphi/PI;
    }

    // First, we calculate Gsigma
    pphi = 0;
    // ... the (1,1) component with cos^2(phi) as prefactor
    sig = 0;
    Gsigma[0] = cquad( wrapphi, params, 0., PI*0.5, relerr*recerr, abserr*recerr);
    Gsigma[0] += cquad( wrapphi, params, PI*0.5, PI, relerr*recerr, abserr*recerr);
    // ... the (2,2) component with sin^2(phi) as prefactor
    sig = 1;
    Gsigma[1] = cquad( wrapphi, params, 0., PI*0.5, relerr*recerr, abserr*recerr);
    Gsigma[1] += cquad( wrapphi, params, PI*0.5, PI, relerr*recerr, abserr*recerr);
    // and the (2,2) component as sin^2 = 1 - cos^2
    sig = 2;
    Gsigma[2] = cquad( wrapphi, params, 0., PI*0.5, relerr*recerr, abserr*recerr);
    Gsigma[2] += cquad( wrapphi, params, PI*0.5, PI, relerr*recerr, abserr*recerr);

    // Second, we calculate Gphi
    pphi = 1;
    Gphi = cquad( wrapphi, params, 0., PI*0.5, relerr*recerr, abserr*recerr);
    Gphi += cquad( wrapphi, params, PI*0.5, PI, relerr*recerr, abserr*recerr);

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
void alpha(double complex alp[Ndim][Ndim], double w, void * p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double wa = (params->wa);
    const double a0 = (params->a0);
    const int muquest = (params->muquest);

    double complex alpinv00,alpinv11,alpinv22,alpinv20, det;
    double complex GI[3][3], GR[3][3];

    alpinv00 =  wa*wa - w*w ;
    if (muquest == 1) {
        alpinv00 = alpinv00 - w*I*mu(w, params);
    }
    alpinv00 = alpinv00/(a0*wa*wa);

    alpinv11 = alpinv00;
    alpinv22 = alpinv00;

    Gint(GI, w, params, 1, 0, 0, 0);
    Gint(GR, w, params, 0, 0, 0, 0);

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

double complex mu(double w, void *p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double gamMu = (params->gamMu);
    
    double complex muresC;
    muresC = gamMu +I*0E0;

    return muresC;
}	

// quantum friction integrand
double IntQF(double w, void * p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const int transroll = (params->transroll);

    double complex GIth[3][3], GIk[3][3], GIkth[3][3];
    double complex alp[3][3], alpI[3][3], S[3][3],  temp1[3][3], temp2[3][3];
    double complex alpdag[3][3];

    /* Creating all needed matrices */
    alpha(alp, w, params);
    dagger(alp, alpdag);
    fancy(alp, alpI, -1);
    Gint(GIth , w, params, 1, 0, 1, 1);
    Gint(GIk  , w, params, 1, 1, 0, 0);
    Gint(GIkth, w, params, 1, 1, 1, 1);
    /* Building the power spectrum S */
    multiply(alp,GIth,temp1);
    multiply(temp1,alpdag,S);

    /* Calculating the trace and return */
    /* We split the calcualtion in translational and rolling friction */
    switch(transroll) {
        case 0:
            GIk[2][0] = 0.;
            GIk[0][2] = 0.;
            GIkth[2][0] = 0.;
            GIkth[0][2] = 0.;
            break;

        case 1:
            GIk[0][0] = 0.;
            GIk[1][1] = 0.;
            GIk[2][2] = 0.;
            GIkth[0][0] = 0.;
            GIkth[1][1] = 0.;
            GIkth[2][2] = 0.;
            break;
    }

    multiply(S,GIk,temp1);
    multiply(alpI,GIkth,temp2);
    /* printf("\nIntQF=%.5e, w=%.5e",(2E0/PI)*creal( -tr(temp1)  + tr(temp2) ),w);*/
    return (2E0/PI)*creal( -tr(temp1)  + tr(temp2) );
}

double QF(double IntQF(), void * p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double abserr = (params->abserr);
    const double relerr = (params->relerr);
    const double wa = (params->wa);
    const double wsp1 = (params->wsp1);
    
    const double v = (params->v);
    const double za = (params->za);
    const double kcut = (params->kcut);
    
    double wcut = kcut*v/(2*za);
    double result, abserr2;
    
    if (wcut < wa*0.9E0) {
	    result = cquad(IntQF, params, 0E0, wcut, relerr, abserr);
             abserr2 = fabs(result)*1E-2;
	    result += cquad(IntQF, params, wcut, 0.9E0*wa, relerr, abserr);
    }
    else {
	    result = cquad(IntQF, params, 0E0, 0.9E0*wa, relerr, abserr);
    }

    abserr2 = fabs(result)*1E-2;
    result += cquad(IntQF, params, 0.9E0*wa, wa, relerr, abserr2);
    
    abserr2 = fabs(result)*1E-2;
    result += cquad(IntQF, params, wa, wsp1, relerr, abserr2);

    abserr2 = fabs(result)*1E-2;
    result += integinf(IntQF, params, wsp1, relerr, abserr2);

    return result;
}

/* analytics */
double F0(void * p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double wsp1 = (params->wsp1);
    const double a0 = (params->a0);
    const double eps0 = (params->eps0);

    return -3*pow(wsp1,5)*a0/(2*PI*eps0);
}

double Fanat(void * p) {

    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double a0 = (params->a0);
    const double g1 = (params->g1);
    const double eps0 = (params->eps0);
    const double wp1 = (params->wp1);
    const double v = (params->v);
    const double za = (params->za);
    const int muquest = (params->muquest);
    const double wa = (params->wa);
    const double beta = (params->beta);

    double result;
    
    result = -63*a0*a0*pow(g1/(eps0*wp1*wp1),2)*pow(v,3)/(pow(PI,3)*pow(2*za,10));

    switch (muquest) {
        case 0:
            break;
        case 1:
            result += -45*creal(mu(0E0, params))*a0*g1*pow(v,3)*pow(2*za,-7)*4*PI/(pow(PI*wa*wp1,2));
            result += -8*v*a0*g1*creal(mu(0E0, params))*4*PI/(pow(wp1*beta*wa,2)*pow(2*za,5));
            break;
    }

    result -= a0*a0*pow(g1/(eps0*wp1*wp1),2)*6*v/(PI*beta*beta*pow(2*za,8));

    return result;
}

double Fanar(void * p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double a0 = (params->a0);
    const double g1 = (params->g1);
    const double eps0 = (params->eps0);
    const double wp1 = (params->wp1);
    const double v = (params->v);
    const double za = (params->za);
    const double beta = (params->beta);

    double result;
    result = 45*a0*a0*pow(g1/(eps0*wp1*wp1),2)*pow(v,3)/(pow(PI,3)*pow(2*za,10));
    result += 3*a0*a0*pow(g1/(eps0*wp1*wp1),2)*v/(PI*beta*beta*pow(2*za,8));
    
    return result;
}

double Ffreet(void * p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;
    const double a0 = (params->a0);
    const double g1 = (params->g1);
    const double eps0 = (params->eps0);
    const double wp1 = (params->wp1);
    const double v = (params->v);
    const double za = (params->za);
    const double beta = (params->beta);

    double result = - a0*a0*pow(g1*4*PI/(eps0*wp1*wp1),2)*3*v/(PI*beta*beta*pow(2*za,8));
    return result;
}

double Ffreer(void * p) {
    /* give parameters needed */
    struct parameters * params = (struct parameters *)p;

    double result = -1./2.*Ffreet(params);
    return result;
}
