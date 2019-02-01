#include "h/plate.h"

int main (int argc, char *argv[]) {
    /* Dummies */
    register unsigned int l; // loop runner
    double spacing, step; // plot parameters
    double QFt, QFr, F0val, Fanarval, Fanatval, Ffreetval, Ffreerval; // dummies

    /* Greetings */
    printf("===========================================\n");
    printf("WELCOME TO QFNUM!\n");
    printf("===========================================\n");

    /* check input file */
    if (argc == 1) {
        printf("No file passed!\n");
        exit(0);
    } else {
        input(argv[1], 1);
    }

    /* create output file */
    // name
    char outfile[strlen(argv[1])+1];
    memset(outfile, '\0', sizeof(outfile)); //clear memory location
    strncpy(outfile, argv[1], strlen(argv[1])-3);
    strcat(outfile, ".out");

    // create file 
    FILE *fp;
    fp = fopen(outfile, "w");
    if (fp == NULL) {
        printf("Couldn't open output file \" %s \" for writing!\n", outfile);
        exit(0);
    };

    /* calculate spacing */
    if (strcmp(inputparams.scale, "log") == 0) {
        spacing = pow(inputparams.stop/inputparams.start,1./inputparams.steps);
    } else {
        printf("Enter valid scaling! (log)\n");
        exit(0);
    };

    /* Starting calculations */
    printf("\n===========================================\n");
    printf("CALCULATION STARTED!\n");

    clock_t c0 = clock();
    clock_t cl0, cl1;
    for (l=0; l<=inputparams.steps; ++l){

        // calculate step
        step = inputparams.start*pow(spacing,l);

        // update running variable
        if (strcmp(inputparams.runvar, "v") == 0) {
            inputparams.v = step;
        } else if (strcmp(inputparams.runvar, "za") == 0) {
            inputparams.za = step/(1.9732705e-7); 
        } else if (strcmp(inputparams.runvar, "T") == 0) {
            inputparams.beta = 1E0/(step/1.16E4);
        } else {
            printf("Enter valid running variable! (v)\n");
            exit(0);
        }

        /* Performing calculations */

        // translational contribution
        inputparams.transroll = 0;

        cl0 = clock();
        QFt = QF(IntQF, &inputparams);
        cl1 = clock();
        printf("Time elapsed QFt: %3.2f sec, ", (double)(cl1-cl0)/CLOCKS_PER_SEC);

        // rolling contribution
        inputparams.transroll = 1;

        cl0 = clock();
        QFr = QF(IntQF, &inputparams);
        cl1 = clock();
        printf("Time elapsed QFr: %3.2f sec\n", (double)(cl1-cl0)/CLOCKS_PER_SEC);

        // analytical limits
        F0val = F0(&inputparams);
        Fanatval = Fanat(&inputparams);
        Fanarval = Fanar(&inputparams);
        Ffreetval = Ffreet(&inputparams);
        Ffreerval = Ffreer(&inputparams);

        /* print results to the screen */
        printf("v        | QFt/F0   | QFr/F0    | Fanat/F0 | Fanar/F0  | Ffreet/F0| Ffreer/F0\n");
        printf("%.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e\n\n",
                inputparams.v, QFt/F0val, QFr/F0val,Fanatval/F0val, Fanarval/F0val,
                Ffreetval/F0val, Ffreerval/F0val);

        /* write results to the file */
        fprintf(fp, "%.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e\n",
                inputparams.v, QFt/F0val, QFr/F0val,Fanatval/F0val, Fanarval/F0val,
                Ffreetval/F0val, Ffreerval/F0val);
        fflush(fp);
    };

    clock_t c1 = clock();
    printf("\n------------------------------\n");

    /* Bye! */
    printf("\n");
    printf("Finished calculating %d points in %3.2f sec. \n", inputparams.steps, (c1-c0)/1.e6);
    printf("\nBYE!\n");
    printf("===========================================\n");

    return 0;
};
