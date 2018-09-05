#include "cyl.h"

int main(void) {
    
    
    /*
     * variables
     */

    double complex omega =0.001;
    omega_p = 1.;
    gamma_p = 0.01;
    c = 1.;
    double h = 0.000001;
    R = 1.;
    eps0 = 1.;
    //double complex x = 2+4*I;
    
   
    /*
     * reflection coefficients
     */

    // calculate  
    double complex result = refCoeffn(1, 1, omega, h);
    
    // print
    printf("Reflection coefficient:\n");
    printf("------------------------\n");
    printf("Ergebnis n: %.18f + %.18f*i \n", creal(result), cimag(result));
    printf("\n");

    /*
     * Green's tensor
     */
    
    // initialize and calculate
    double complex g[3][3];
    int i,j;
    greencent(omega, h, g); 

    // print tensor
    printf("Central Green's tensor:\n");
    printf("-----------------------\n");
    for (i = 0; i <= 2; i++)
    {
        for (j = 0; j <= 2; j++)
           printf("%.8f + %.8f*i   ", creal(g[i][j]), cimag(g[i][j]));
        printf("\n\n");
    } 

    double complex rho = 0.001;
    int N = 1;
    double complex gtest[3][3];
    greenfull(N, omega, h, rho, gtest); 
    
    // print tensor
    printf("Full Green's tensor:\n");
    printf("--------------------\n");
    for (i = 0; i <= 2; i++)
    {
        for (j = 0; j <= 2; j++)
           printf("%.8f + %.8f*i   ", creal(gtest[i][j]), cimag(gtest[i][j]));
        printf("\n\n");
    } 
};
