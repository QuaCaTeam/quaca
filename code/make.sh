gcc -Wall -g -I/opt/vc/include -c QF.c
gcc -g -Wall -L/opt/vc/lib QF.o -lgsl -lgslcblas -lm -o qf

