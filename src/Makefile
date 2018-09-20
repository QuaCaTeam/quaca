#
# Makefile 
#

# directories
SRCDIR = ../src
INCDIR = ../src/h

# options
CC = gcc
CFLAGS = -Wall # show warnings 
INCLUDES = -I/usr/include -I$(INCDIR) # include headers
LFLAGS = -L/usr/lib 
LIBS = -lgsl -lgslcblas -lm -lflint-arb -lflint -lmpfr # libraries (gsl and arb), on arch -larb on debian -lflint-arb

# file
SRC = qfPlate.c cyl.c plate.c qfhelp.c
OBJS = ${SRC:.c=.o}
PROG = qfPlate

all: $(PROG) clean

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# make docs
docs:
	@doxygen doc/doxconf

clean:
	rm *.o *.bak

.PHONY: all clean tester
