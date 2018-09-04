#include "green.h"

int main(void) {
    
    
    /*
     * variables
     */

    double complex omega =0.001;
    double omega_p = 1;
    double gamma = 0.01;
    double c = 1;
    double h = 0.000001;
    double R = 1;
    //double complex x = 2+4*I;
    
   
    /*
     * reflection coefficients
     */

    // calculate  
    double complex result = rMM1(omega, omega_p, gamma, c, h, R);
    double complex result2 = refCoeffn(1, 1, omega, omega_p, gamma, c, h, R);
    
    // print
    printf("Reflection coefficient:\n");
    printf("------------------------\n");
    printf("Ergebnis 1: %.18f + %.18f*i \n", creal(result), cimag(result));
    printf("Ergebnis n: %.18f + %.18f*i \n", creal(result2), cimag(result2));
    printf("\n");

    /*
     * Green's tensor
     */
    
    // initialize and calculate
    double complex g[2][2];
    int i,j;
    greencent(omega, omega_p, gamma, c, h, R, 1, g); 

    // print tensor
    printf("Green's tensor:\n");
    printf("---------------\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
           printf("%f + %f*i   ", creal(g[i][j]), cimag(g[i][j]));
        printf("\n\n");
    } 
};
