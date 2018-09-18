#
# Makefile 
#

# options
CC = gcc
CFLAGS = -Wall -I/usr/include # show warnings and include
LDFLAGS = -L/usr/lib -lgsl -lgslcblas -lm -lflint-arb -lflint -lmpfr # libraries (gsl and arb), on arch -larb on debian -lflint-arb

# file
SRC = tester.c ref.c QF.c Omega.c
OBJ = ${SRC:.c=.o}
PROG = ${SRC:.c=}

all: $(PROG) clean

# make object files
.c.o:
	$(CC) -c $(CFLAGS) $<

# link object files, create executable
$(PROG): %: %.o
	$(CC) -o $@ $< $(LDFLAGS)

# make docs
docs:
	@cd ../docs && doxygen doxconf

clean:
	rm *.o *.bak

.PHONY: all clean tester
