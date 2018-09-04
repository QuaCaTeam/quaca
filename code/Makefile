#
# Makefile for testing
#

# options
CC = gcc
CFLAGS = -Wall -I/usr/include # show warnings and include
LDFLAGS = -L/usr/lib -lgsl -lgslcblas -lm -lflint-arb -lflint -lmpfr # libraries (gsl and arb), on arch -larb on debian -lflint-arb

# file
SRC = tester.c
OBJ = ${SRC:.c=.o}
PROG = tester

all: $(PROG) clean

# make object files
.c.o:
	$(CC) -c $(CFLAGS) $<

# link object files, create executable
$(PROG): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o

.PHONY: all clean tester
