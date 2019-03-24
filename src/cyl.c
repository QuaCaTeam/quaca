/* --------- */
/* LIBRARIES */
/* --------- */

#include "h/qfhelp.h"
#include "h/cyl.h"

/* --------- */
/* FUNCTIONS */
/* --------- */

// input routine
void inputCyl(char file[], unsigned int verbose) {
    FILE * fr = fopen(file, "rt");

    if (fr == NULL) {
        printf("file %s not found", file);
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
                    inputparamsCyl.wp1 = i;
                } else if (strcmp(s,"g1")==0) {
                    inputparamsCyl.g1 = i;
                } else if (strcmp(s,"einf")==0) {
                    inputparamsCyl.einf = i;
                } else if (strcmp(s,"v")==0) {
                    inputparamsCyl.v = i;
                } else if (strcmp(s,"R")==0) {
                    inputparamsCyl.R = i;
                } else if (strcmp(s,"T")==0) {
                    inputparamsCyl.beta = 1E0/(i*kB);
                } else if (strcmp(s,"a0")==0) {
                    inputparamsCyl.a0 = i;
                } else if (strcmp(s,"kcut")==0) {
                    inputparamsCyl.kcut = i;
                } else if (strcmp(s,"relerr")==0) {
                    inputparamsCyl.relerr = i;
                } else if (strcmp(s,"abserr")==0) {
                    inputparamsCyl.abserr = i;
                } else if (strcmp(s,"recerr")==0) {
                    inputparamsCyl.recerr = i;
                } else if (strcmp(s,"wa")==0) {
                    inputparamsCyl.wa = i;
                } else if (strcmp(s,"vF")==0) {
                    inputparamsCyl.vF = i;
                } else if (strcmp(s,"aF")==0) {
                    inputparamsCyl.aF = i;
                } else if (strcmp(s,"gamMu")==0) {
                    inputparamsCyl.gamMu = i;
                } else if (strcmp(s,"muquest")==0) {
                    inputparamsCyl.muquest = (int) i;
                } else if (strcmp(s,"start")==0) {
                    inputparamsCyl.start = i;
                } else if (strcmp(s,"stop")==0) {
                    inputparamsCyl.stop = i;
                } else if (strcmp(s,"steps")==0) {
                    inputparamsCyl.steps = (int) i;
                } else {
                    printf("Unrecongized parameter : \"%s\"\n", s);
                }
            }

            if ( sscanf(line, "runvar %s", s) == 1 ) {
                strcpy(inputparamsCyl.runvar, s); 
            } 

            if ( sscanf(line, "scale %s", s) == 1 ) {
                strcpy(inputparamsCyl.scale, s); 
            } 
        }
    }


    // recall the read parameters
    if (verbose == 1) {
        printf("INPUT PARAMETERS:\n");
        printf("// atomic constants\n");
        printf("v   = %.5e\n",inputparamsCyl.v);
        printf("R  = %.5e\n",inputparamsCyl.R);
        printf("wa  = %.5e\n",inputparamsCyl.wa);
        printf("a0  = %.5e\n",inputparamsCyl.a0);
        printf("\n// material parameters\n");
        printf("einf = %.5e\n",inputparamsCyl.einf);
        printf("wp1  = %.5e\n",inputparamsCyl.wp1);
        printf("g1   = %.5e\n",inputparamsCyl.g1);
        printf("vF   = %.5e\n",inputparamsCyl.vF);
        printf("aF   = %.5e\n",inputparamsCyl.aF);
        printf("T    = %.5e\n",1.16e4/inputparamsCyl.beta);
        printf("\n// numerical specifications\n");
        printf("kcut   = %.5e\n",inputparamsCyl.kcut);
        printf("relerr = %.5e\n",inputparamsCyl.relerr);
        printf("recerr = %.5e\n",inputparamsCyl.recerr);
        printf("abserr = %.5e\n",inputparamsCyl.abserr);
        printf("\n\nCalculating %s from %.5e to %.5e with %d points on %s scale.\n", inputparamsCyl.runvar, inputparamsCyl.start, inputparamsCyl.stop, inputparamsCyl.steps, inputparamsCyl.scale);
    }

    // transforming the remaining SI units to natural units
    // or other convenient forms
    //inputparamsCyl.R   = inputparamsCyl.R/(hbar*cvac);
    inputparamsCyl.wsp1 = inputparamsCyl.wp1/sqrt(1.+inputparamsCyl.einf);
    inputparamsCyl.eps0= 1./(4*PI);
}

const double complex hankel1(double complex n, double complex x) {
    double complex result = ac_besselj(n,x) + I*ac_bessely(n,x);
    return result;
};

const double complex hn(int n, double complex x) {
    double complex result = (n/x*hankel1(n,x) - hankel1(n+1,x))/hankel1(n,x);
    return result;
};

const double complex besseljn(int n, double complex x) {
    double complex result = (n/x*ac_besselj(n,x) - ac_besselj(n+1,x))/ac_besselj(n,x);
    return result;
};

double complex rNN(int n, double w, double h, void * p) {
    /* give parameters needed */
    struct parametersCyl * params = (struct parametersCyl *)p;
    const double einf = (params->einf);
    const double wp1 = (params->wp1);
    const double g1 = (params->g1);
    const double R = (params->R);

    double complex epsw, x, x1;
    double complex A, BNN, C;
    double complex result;

    epsw = einf - wp1*wp1/(w*w + w*g1*I);
    x = R*csqrt(w*w - h*h);
    x1 = csqrt(w*w*epsw - h*h);

    A = - n*n*w*w*h*h*pow(R,4)*(epsw - 1)*(epsw-1);
    BNN = x1*x1*x*x*(epsw*hn(n,x1)*hn(n,x1)*x*x - x*x1*(epsw*hn(n,x1)*besseljn(n,x) + hn(n,x1)*hn(n,x)) + hn(n,x)*besseljn(n,x)*x1*x1);
    C = x1*x1*x*x*(epsw*hn(n,x1)*hn(n,x1)*x*x - (epsw + 1)*hn(n,x1)*besseljn(n,x)*x1*x + besseljn(n,x)*besseljn(n,x)*x1*x1);
    
    result = - hankel1(n,x)/ac_besselj(n,x)*(A + BNN)/(A + C);
    
    return result;
};

double complex rMM(int n, double w, double h, void * p) {
    /* give parameters needed */
    struct parametersCyl * params = (struct parametersCyl *)p;
    const double einf = (params->einf);
    const double wp1 = (params->wp1);
    const double g1 = (params->g1);
    const double R = (params->R);

    double complex epsw, x, x1;
    double complex A, BMM, C;
    double complex result;

    epsw = einf - wp1*wp1/(w*w + w*g1*I);
    x = R*csqrt(w*w - h*h);
    x1 = csqrt(w*w*epsw - h*h);

    A = - n*n*w*w*h*h*pow(R,4)*(epsw - 1)*(epsw-1);
    BMM = x1*x1*x*x*(epsw*hn(n,x1)*hn(n,x1)*x*x - x*x1*(hn(n,x1)*besseljn(n,x) + epsw*hn(n,x1)*hn(n,x)) + hn(n,x)*besseljn(n,x)*x1*x1);
    C = x1*x1*x*x*(epsw*hn(n,x1)*hn(n,x1)*x*x - (epsw + 1)*hn(n,x1)*besseljn(n,x)*x1*x + besseljn(n,x)*besseljn(n,x)*x1*x1);
    
    result = - hankel1(n,x)/ac_besselj(n,x)*(A + BMM)/(A + C);
    
    return result;
};


double complex rMN(int n, double w, double h, void * p) {
    /* give parameters needed */
    struct parametersCyl * params = (struct parametersCyl *)p;
    const double einf = (params->einf);
    const double wp1 = (params->wp1);
    const double g1 = (params->g1);
    const double R = (params->R);

    double complex epsw, x, x1;
    double complex A, BMN, C;
    double complex result;

    epsw = einf - wp1*wp1/(w*w + w*g1*I);
    x = R*csqrt(w*w - h*h);
    x1 = csqrt(w*w*epsw - h*h);

    A = - n*n*w*w*h*h*pow(R,4)*(epsw - 1)*(epsw-1);
    BMN = I*n*x1*x1*x*R*w*h*R*(1 - epsw)*(hn(n,x) - besseljn(n,x));
    C = x1*x1*x*x*(epsw*hn(n,x1)*hn(n,x1)*x*x - (epsw + 1)*hn(n,x1)*besseljn(n,x)*x1*x + besseljn(n,x)*besseljn(n,x)*x1*x1);
    
    result = - hankel1(n,x)/ac_besselj(n,x)*(BMN)/(A + C);
    
    return result;
};


double complex rNNNF(int n, double w, double h, void * p) {
    /* give parameters needed */
    struct parametersCyl * params = (struct parametersCyl *)p;
    const double einf = (params->einf);
    const double wp1 = (params->wp1);
    const double g1 = (params->g1);
    const double R = (params->R);

    double complex epsw, result, fn, x;

    x = h*R;

    epsw = einf - wp1*wp1/(w*w + w*g1*I);

    fn = - (ac_besselk(n,x)*0.5*(ac_besseli(n-1,x)+ac_besseli(n+1,x)))/(ac_besseli(n,x)*(-0.5)*(ac_besselk(n-1,x) + ac_besselk(n+1,x)));

    result = pow(-1,n)*2*I/PI*ac_besselk(n,x)/ac_besseli(n,x)*(epsw -1)/(epsw + fn);

    return result;
};

void GCNFint(double Gten[Ndim][Ndim], double w, void * p, int RorI) {
    /* give parameters needed */
    struct parametersCyl * params = (struct parametersCyl *)p;
    const double v = (params->v);
    const double eps0 = (params->eps0);
    const double relerr = (params->relerr);
    const double kcut = (params->kcut);
    const double abserr = (params->abserr);

    double Gten1, Gten2;
    register unsigned int comp;
    
    double wraph(double h) {
       double complex res = 0;
       double resh = 0;
       double wdopl = w+v*h;

       switch(comp) {
           case 0:
               res = I*h*h/(8*eps0)*rNNNF(1,wdopl,h,params);
               break;
           case 1:
               res = -2*I*h*h/(8*eps0)*rNNNF(0,wdopl,h,params);
               break;
       }

        
       switch(RorI) {
           case 0:
               resh = creal(res);
               break;
           case 1:
               resh = cimag(res);
               break;
       }

       return resh;
    };

    comp = 0;
    Gten1 = cquad(wraph, params, 0, kcut, relerr, abserr);
    
    comp = 1;
    Gten2 = cquad(wraph, params, 0, kcut, relerr, abserr);

    Gten[0][0] = Gten1;
    Gten[1][1] = Gten1;
    Gten[2][2] = Gten2;
    Gten[0][2] = 0.;
    Gten[2][0] = 0.;
    Gten[0][1] = 0.;
    Gten[1][0] = 0.;
    Gten[2][1] = 0.;
    Gten[1][2] = 0.;
};


void GCint(double Gten[Ndim][Ndim], double w, void * p, int RorI) {
    /* give parameters needed */
    struct parametersCyl * params = (struct parametersCyl *)p;
    const double v = (params->v);
    const double eps0 = (params->eps0);
    const double relerr = (params->relerr);
    const double abserr = (params->abserr);
    const double kcut = (params->kcut);

    double Gten1, Gten2;
    register unsigned int comp;
    
    double wraph(double h) {
       double complex res = 0;
       double resh = 0;
       double wdopl = w+v*h;

       switch(comp) {
           case 0:
               res = I/(8*eps0)*(wdopl*wdopl*rMM(1,wdopl,h,params) + h*h*rNN(1,wdopl,h,params) - 2*I*wdopl*h*rMN(1,wdopl,h,params));
               break;
           case 1:
               res = I/(4*eps0)*csqrt(wdopl*wdopl - h*h)*rNN(0,wdopl,h,params);
               break;
       }

        
       switch(RorI) {
           case 0:
               resh = creal(res);
               break;
           case 1:
               resh = cimag(res);
               break;
       }

       return resh;
    };

    comp = 0;
    Gten1 = cquad(wraph, params, 0, kcut, relerr, abserr);
    
    comp = 1;
    Gten2 = cquad(wraph, params, 0, kcut, relerr, abserr);

    Gten[0][0] = Gten1;
    Gten[1][1] = Gten1;
    Gten[2][2] = Gten2;
    Gten[0][2] = 0.;
    Gten[2][0] = 0.;
    Gten[0][1] = 0.;
    Gten[1][0] = 0.;
    Gten[2][1] = 0.;
    Gten[1][2] = 0.;
};
