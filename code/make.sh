gcc -Wall -g -I/opt/vc/include -c QF.c
gcc -g -Wall -L/opt/vc/lib QF.o -lgsl -lgslcblas -lm -o qf -O2

gcc -Wall -g -I/opt/vc/include -c Omega.c
gcc -g -Wall -L/opt/vc/lib Omega.o -lgsl -lgslcblas -lm -o om -O2

gcc -Wall -g -I/opt/vc/include -c IntOmega.c
gcc -g -Wall -L/opt/vc/lib IntOmega.o -lgsl -lgslcblas -lm -o intom -O2

gcc -Wall -g -I/opt/vc/include -c ref.c
gcc -g -Wall -L/opt/vc/lib ref.o -lgsl -lgslcblas -lm -o ref -O2

