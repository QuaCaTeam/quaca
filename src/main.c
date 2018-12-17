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
// Main program

int main (int argc, char *argv[]) {

    /* Plot parameters */
    int maxi;                        // plot points
    double sta;                      // start value of the calculation,
    double sto;                      // final value of the calculation
    double spac;                     // and the respective spacing

    /* Flag */

    /* Dummies */
    int l;
    double QFt, QFr;// QFtfree, QFrfree;
    double F0, Fanar, Fanat, Ffreet, Ffreer;
    double eps0= 1./(4*PI);

    // greetings
    printf("===========================================\n");
    printf("WELCOME TO QFNUM!\n");
    printf("===========================================\n");

    // input parameters
    input(1);    

    sta = 1E-3; 
    //printf("Start velocity [in c]: ");
    //scanf("%.10e",&sta);  
    //printf("Read sta = %.10e",sta);

    sto = 1E-2;
    //printf("Stop velocity [in c]: ");
    //scanf("%.10e",&sto);  
    //printf("Read sto = %.10e",sto);
 
    maxi = 100;
    //printf("Number of point: ");
    //scanf("%d",&maxi);  
    //printf("Read maxi = %d",maxi);

    spac = pow(sto/sta,1./maxi);

    /* print evaluation regime */
    printf("v= %.5e to %.5e\n", sta, sto);
    printf("with %5i points\n", maxi);

    /* Starting calculations */
    printf("===========================================\n");
    printf("CALCULATION STARTED!\n");
    printf("v               | QFt/F0          | QFr/F0           | Fanat/F0        | Fanar/F0         | Ffreet/F0       | Ffreer/F0\n");

    clock_t c0 = clock();
    for (l=0; l<=maxi; ++l){
        v = sta*pow(spac,l);

    /* Performing calculations */
	 /* translational contribution */
 	 transroll = 0;
  		QFt =  integ(IntQF,0E0,0.9E0*wa,relerr, abserr);
  		abserr = fabs(QFt)*1E-2;
  		QFt += integ(IntQF,0.9E0*wa,wa,relerr, abserr);
  		abserr = fabs(QFt)*1E-2;
  		QFt += integ(IntQF,wa,wsp1,relerr, abserr);
		abserr = fabs(QFt)*1E-2;
 		QFt += integinf(IntQF,wsp1,relerr, abserr);
  		abserr = 1E-200;
  	/* rolling contribution */
  	transroll = 1;
  		QFr =  integ(IntQF,0E0,0.9E0*wa,relerr, abserr);
  		abserr = fabs(QFr)*1E-2;
 		QFr += integ(IntQF,0.9E0*wa,wa,relerr, abserr);
		abserr = fabs(QFr)*1E-2;
		QFr += integ(IntQF,wa,wsp1,relerr, abserr);
		abserr = fabs(QFr)*1E-2;
		QFr += integinf(IntQF,wsp1,relerr, abserr);
		abserr = 1E-200;
  /* Calculate normalization constant */
  	F0 = -3*pow(wsp1,5)*a0/(2*PI*eps0);
  /* Calculating analytical approximation for small velocities */
  	Fanat = -63*a0*a0*pow(g1/(eps0*wp1*wp1),2)*pow(v,3)/(pow(PI,3)*pow(2*za,10));
  	Fanar = -45*Fanat/63.;
  	if (muquest == 1){
   	 Fanat += -45*creal(mu(0E0))*a0*g1*pow(v,3)*pow(2*za,-7)*4*PI/(pow(PI*wa*wp1,2));
   	 Fanat += -8*v*a0*g1*creal(mu(0E0))*4*PI/(pow(wp1*beta*wa,2)*pow(2*za,5));
  	}
  	Fanat = Fanat - a0*a0*pow(g1/(eps0*wp1*wp1),2)*6*v/(PI*beta*beta*pow(2*za,8));
  	Fanar = Fanar + 3*a0*a0*pow(g1/(eps0*wp1*wp1),2)*v/(PI*beta*beta*pow(2*za,8));
  	Ffreet= - a0*a0*pow(g1*4*PI/(eps0*wp1*wp1),2)*3*v/(PI*beta*beta*pow(2*za,8));
  	Ffreer= -Ffreet/2.;
  
  /* Print result to the screen */
   printf("\n");
   printf("v= %.5e\n", v);
   printf("F= %.5e\n", (QFt+QFr)/F0 );
  /* write to the file */
   fprintf(fp, "%.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e\n",
             v, QFt/F0, QFr/F0,Fanat/F0, Fanar/F0,
             Ffreet/F0, Ffreer/F0);
   fflush(fp);
  /* Progress bar */
   prog = l/(double)maxi;
   printProg(prog);
  }
clock_t c1 = clock();
printf("\n------------------------------\n");

/* Bye! */
printf("\n");
printf("Finished calculating %d points in %3.2f sec. \n", maxi, (c1-c0)/1.e6);
printf("\nBYE!\n");
printf("===========================================\n");

return 0;
};

////////////////////////////////////////////////////////////////////////////////
