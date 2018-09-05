#include "cyl.h"

int main(void) {
    
    
    /*
     * variables
     */

    double complex omega = 0.001;
    double omega_p = 1;
    double gamma = 9;
    double c = 1;
    double h = 0.000001;
    double R = 1;
    double eps0 = 1;
    //double complex x = 2+4*I;
    
   
    /*
     * reflection coefficients
     */

    // calculate  
    double complex result = refCoeffn(1, 1, omega, omega_p, gamma, c, h, R);
    
    // print
    printf("Reflection coefficient:\n");
    printf("------------------------\n");
    printf("Ergebnis n: %.18f + %.18f*i \n", creal(result), cimag(result));
    printf("\n");

    /*
     * Green's tensor
     */
    
    // initialize and calculate
    double complex g[2][2];
    int i,j;
    greencent(omega, omega_p, gamma, c, h, R, eps0, g); 

    // print tensor
    printf("Central Green's tensor:\n");
    printf("-----------------------\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
           printf("%.8f + %.8f*i   ", creal(g[i][j]), cimag(g[i][j]));
        printf("\n\n");
    } 

    double complex rho = 0.001;
    int N = 0;
    double complex gtest[2][2];
    greenfull(N, omega, omega_p, gamma, c, h, R, eps0, rho, gtest); 
    
    // print tensor
    printf("Full Green's tensor:\n");
    printf("--------------------\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
           printf("%.8f + %.8f*i   ", creal(gtest[i][j]), cimag(gtest[i][j]));
        printf("\n\n");
    } 
};
