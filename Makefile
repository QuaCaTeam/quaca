#
# Makefile 
#

# options
CC = gcc
CFLAGS = -Wall -I/usr/include # show warnings and include
LIBS = -L/usr/lib -lgsl -lgslcblas -lm -lflint-arb -lflint -lmpfr # libraries (gsl and arb), on arch -larb on debian -lflint-arb

# directories
SRCDIR = src
HEADDIR = /header

# file
SRC = tester.c ref.c QF.c Omega.c IntOmega.c
OBJ = ${SRC:.c=.o}
PROG = ${SRC:.c=}

all: $(PROG) clean

# make object files
.c.o:
	$(CC) -c $(CFLAGS) $<

# link object files, create executable
$(PROG): %: %.o
	$(CC) -o $@ $< $(LIBS)

# make docs
docs:
	@doxygen doc/doxconf

clean:
	rm *.o *.bak

.PHONY: all clean tester
